#!/usr/bin/env python3
"""
ConformSeek API Server
======================
FastAPI backend providing:
  1. PROPKA pKa analysis — accepts PDB ID or structure file, returns per-residue pKa values
  2. BRENDA pH optimum lookup — queries the BRENDA SOAP API (Bring Your Own Key)
     for experimental pH optima for a given EC number

Run with:
    pip install fastapi uvicorn propka biopython numpy pandas requests zeep
    uvicorn conformseek_server:app --reload --port 8000

BRENDA Setup (one-time):
  1. Register free at https://www.brenda-enzymes.org/register.php
  2. POST /brenda/credentials with {"email": "...", "password": "..."}
  3. Query pH data: GET /brenda/3.2.1.17

The React frontend calls these endpoints to replace the random pKa placeholders
with real PROPKA predictions, and to get high-confidence BRENDA pH data.

Author: Calvin (ConformSeek project)
"""

import io
import os
import json
import re
import tempfile
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import requests
from fastapi import FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# —— PROPKA import ——
try:
    import propka.run
    import propka.molecular_container
    HAS_PROPKA = True
except ImportError:
    HAS_PROPKA = False
    print("WARNING: PROPKA not installed. pip install propka")

# —— BioPython import (for mmCIF → PDB conversion) ——
# NOTE: three_to_one moved in BioPython 1.80+; handle both locations gracefully
HAS_BIOPYTHON = False
try:
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO
    from Bio.PDB.Polypeptide import is_aa
    HAS_BIOPYTHON = True
    try:
        from Bio.PDB.Polypeptide import three_to_one
    except ImportError:
        try:
            from Bio.Data.IUPACData import protein_letters_3to1 as _p3to1
            def three_to_one(three):
                return _p3to1.get(three.upper(), "X")
        except ImportError:
            def three_to_one(three):
                _MAP = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
                        "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
                        "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
                        "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
                return _MAP.get(three.upper(), "X")
except ImportError:
    print("WARNING: BioPython not installed. pip install biopython")

# —— BRENDA SOAP API (zeep) ——
try:
    from zeep import Client as ZeepClient, Settings as ZeepSettings
    HAS_ZEEP = True
except ImportError:
    HAS_ZEEP = False
    print("WARNING: zeep not installed. pip install zeep")
    print("         BRENDA SOAP API lookup will be unavailable.")


# ═══════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════

LOG = logging.getLogger("conformseek.api")
logging.basicConfig(level=logging.INFO)

# —— BRENDA SOAP API configuration ——
BRENDA_WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

# Credentials storage directory
CONFORMSEEK_CONFIG_DIR = Path.home() / ".conformseek"
CONFORMSEEK_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
BRENDA_CREDS_FILE = CONFORMSEEK_CONFIG_DIR / "brenda_credentials.json"

# In-memory SOAP client + credential cache
_brenda_soap_client: Optional[object] = None
_brenda_credentials: Optional[Dict[str, str]] = None  # {"email": ..., "password_hash": ...}

# In-memory cache: EC number → list of {ph, organism}  (populated per-query, persisted to disk)
BRENDA_QUERY_CACHE_FILE = CONFORMSEEK_CONFIG_DIR / "brenda_soap_cache.json"
_brenda_query_cache: Dict[str, List[dict]] = {}

# Standard textbook pKa values (fallback)
STANDARD_PKA = {
    "ASP": 3.65, "GLU": 4.25, "HIS": 6.00,
    "CYS": 8.18, "LYS": 10.53, "TYR": 10.07,
}
TITRATABLE_TYPE = {
    "ASP": "acid", "GLU": "acid", "HIS": "base",
    "CYS": "acid", "LYS": "base", "TYR": "acid",
}
TITRATABLE_SET = set(STANDARD_PKA.keys())


# ═══════════════════════════════════════════════════════════════
# FastAPI app
# ═══════════════════════════════════════════════════════════════

app = FastAPI(
    title="ConformSeek API",
    description="Backend for PROPKA pKa analysis and BRENDA pH lookup (SOAP API + BYOK)",
    version="2.0.0",
)

# Allow React dev server (localhost:3000) to call this API
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "http://localhost:3001",
        "http://127.0.0.1:3001",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ═══════════════════════════════════════════════════════════════
# Pydantic models
# ═══════════════════════════════════════════════════════════════

class PropkaResidue(BaseModel):
    pos: int
    one: str
    three: str
    chain: str
    pka: float
    type: str          # "acid" or "base"
    cls: str           # "Acidic", "Basic", "Thiol", "Phenolic"
    std_pka: float     # standard textbook pKa for comparison
    pka_shift: float   # environment-induced shift from standard

class PropkaResult(BaseModel):
    pdb_id: str
    chain_id: str
    n_residues: int
    n_titratable: int
    titratable: List[PropkaResidue]
    pI: float
    method: str        # "propka" or "standard_fallback"

class BrendaResult(BaseModel):
    ec_number: str
    ph_values: List[float]
    ph_mean: float
    ph_median: float
    ph_min: float
    ph_max: float
    n_entries: int
    organisms: List[str]

class HealthResponse(BaseModel):
    status: str
    propka_available: bool
    biopython_available: bool
    brenda_soap_available: bool
    brenda_credentials_set: bool
    brenda_cached_ec_count: int


# ═══════════════════════════════════════════════════════════════
# BRENDA SOAP API client + credential management (BYOK)
# ═══════════════════════════════════════════════════════════════

import hashlib


def _load_brenda_credentials() -> Optional[Dict[str, str]]:
    """Load saved BRENDA credentials from disk."""
    global _brenda_credentials
    if _brenda_credentials is not None:
        return _brenda_credentials

    if BRENDA_CREDS_FILE.exists():
        try:
            with open(BRENDA_CREDS_FILE) as f:
                creds = json.load(f)
            if creds.get("email") and creds.get("password_hash"):
                _brenda_credentials = creds
                LOG.info(f"Loaded BRENDA credentials for {creds['email']}")
                return _brenda_credentials
        except Exception as e:
            LOG.warning(f"Failed to load BRENDA credentials: {e}")

    return None


def _save_brenda_credentials(email: str, password: str):
    """Save BRENDA credentials to disk (password stored as SHA-256 hash)."""
    global _brenda_credentials, _brenda_soap_client
    password_hash = hashlib.sha256(password.encode("utf-8")).hexdigest()
    _brenda_credentials = {"email": email, "password_hash": password_hash}
    _brenda_soap_client = None  # reset client so it reconnects

    try:
        with open(BRENDA_CREDS_FILE, "w") as f:
            json.dump(_brenda_credentials, f)
        # Restrict file permissions (owner read/write only)
        os.chmod(BRENDA_CREDS_FILE, 0o600)
        LOG.info(f"Saved BRENDA credentials for {email}")
    except Exception as e:
        LOG.warning(f"Failed to save BRENDA credentials: {e}")


def _get_brenda_soap_client():
    """Get or create a zeep SOAP client for BRENDA."""
    global _brenda_soap_client
    if _brenda_soap_client is not None:
        return _brenda_soap_client

    if not HAS_ZEEP:
        raise HTTPException(503, "zeep library not installed. Run: pip install zeep")

    try:
        settings = ZeepSettings(strict=False, xml_huge_tree=True)
        _brenda_soap_client = ZeepClient(BRENDA_WSDL, settings=settings)
        LOG.info("BRENDA SOAP client initialized")
        return _brenda_soap_client
    except Exception as e:
        LOG.error(f"Failed to initialize BRENDA SOAP client: {e}")
        raise HTTPException(503, f"Failed to connect to BRENDA SOAP service: {e}")


def _query_brenda_ph_optimum(ec_number: str) -> List[dict]:
    """
    Query BRENDA SOAP API for pH optimum data for a given EC number.

    Uses the getPhOptimum method with the user's BYOK credentials.
    Results are cached locally to avoid redundant API calls.

    Returns list of {ph: float, organism: str}
    """
    # Check query cache first
    if ec_number in _brenda_query_cache:
        LOG.info(f"BRENDA cache hit for EC {ec_number}")
        return _brenda_query_cache[ec_number]

    creds = _load_brenda_credentials()
    if creds is None:
        raise HTTPException(
            401,
            "BRENDA credentials not configured. Use POST /brenda/credentials "
            "to set your BRENDA email and password. Register free at "
            "https://www.brenda-enzymes.org/register.php"
        )

    client = _get_brenda_soap_client()

    try:
        # BRENDA SOAP getPhOptimum method signature:
        # (email, password_hash, "ecNumber*X.X.X.X", "phOptimum*",
        #  "phOptimumMaximum*", "commentary*", "organism*", "literature*")
        parameters = (
            creds["email"],
            creds["password_hash"],
            f"ecNumber*{ec_number}",
            "phOptimum*",
            "phOptimumMaximum*",
            "commentary*",
            "organism*",
            "literature*",
        )
        result = client.service.getPhOptimum(*parameters)
    except Exception as e:
        error_str = str(e).lower()
        if "credentials" in error_str or "password" in error_str or "access" in error_str:
            raise HTTPException(
                401,
                "BRENDA authentication failed. Check your email and password "
                "via POST /brenda/credentials. Your password may have changed, "
                "or your account may not be activated yet."
            )
        LOG.error(f"BRENDA SOAP query failed for EC {ec_number}: {e}")
        raise HTTPException(502, f"BRENDA SOAP API error: {e}")

    # Parse the SOAP response string into structured entries
    entries = _parse_soap_ph_response(result, ec_number)

    # Cache the result
    _brenda_query_cache[ec_number] = entries
    _save_query_cache()

    return entries


def _parse_soap_ph_response(result: str, ec_number: str) -> List[dict]:
    """
    Parse the BRENDA SOAP getPhOptimum response.

    BRENDA returns a string with entries separated by '!' and fields
    separated by '#' in key*value format. Example fields:
      ecNumber*1.1.1.1  organism*Homo sapiens  phOptimum*7.5
      commentary*...  literature*...
    """
    entries = []
    if not result:
        return entries

    # Handle both string and list responses
    if isinstance(result, str):
        raw_entries = result.split("!")
    elif isinstance(result, (list, tuple)):
        raw_entries = result
    else:
        raw_entries = [str(result)]

    for raw in raw_entries:
        raw = str(raw).strip()
        if not raw:
            continue

        # Parse key*value pairs separated by #
        fields = {}
        for part in raw.split("#"):
            part = part.strip()
            if "*" in part:
                key, _, value = part.partition("*")
                fields[key.strip()] = value.strip()

        # Extract pH value and organism
        ph_str = fields.get("phOptimum", "")
        organism = fields.get("organism", "unknown")

        if not ph_str:
            continue

        try:
            # Handle range values like "6.5-7.0" → take midpoint
            if "-" in ph_str and not ph_str.startswith("-"):
                parts = ph_str.split("-", 1)
                lo, hi = float(parts[0]), float(parts[1])
                ph_val = (lo + hi) / 2.0
            else:
                ph_val = float(ph_str)

            if 0 <= ph_val <= 14:
                entries.append({
                    "ph": round(ph_val, 2),
                    "organism": organism,
                })
        except (ValueError, TypeError, IndexError):
            continue

    LOG.info(f"BRENDA SOAP: EC {ec_number} → {len(entries)} pH optimum entries")
    return entries


def _load_query_cache():
    """Load the query cache from disk on startup."""
    global _brenda_query_cache
    if BRENDA_QUERY_CACHE_FILE.exists():
        try:
            with open(BRENDA_QUERY_CACHE_FILE) as f:
                _brenda_query_cache = json.load(f)
            LOG.info(f"Loaded BRENDA query cache: {len(_brenda_query_cache)} EC numbers")
        except Exception as e:
            LOG.warning(f"Failed to load BRENDA query cache: {e}")
            _brenda_query_cache = {}


def _save_query_cache():
    """Persist the query cache to disk."""
    try:
        with open(BRENDA_QUERY_CACHE_FILE, "w") as f:
            json.dump(_brenda_query_cache, f)
    except Exception as e:
        LOG.warning(f"Failed to save BRENDA query cache: {e}")


# ═══════════════════════════════════════════════════════════════
# PROPKA runner
# ═══════════════════════════════════════════════════════════════

THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

CLS_MAP = {
    "ASP": "Acidic", "GLU": "Acidic", "HIS": "Basic",
    "CYS": "Thiol", "LYS": "Basic", "TYR": "Phenolic",
}


def _download_pdb_file(pdb_id: str) -> str:
    """Download a PDB file from RCSB and return path to temp file."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        # Try mmCIF and convert
        url_cif = f"https://files.rcsb.org/download/{pdb_id.lower()}.cif"
        resp = requests.get(url_cif, timeout=30)
        if resp.status_code != 200:
            raise HTTPException(404, f"Could not download structure for {pdb_id}")

        # Convert mmCIF → PDB
        if not HAS_BIOPYTHON:
            raise HTTPException(500, "BioPython required for mmCIF conversion")

        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False, mode="w") as f:
            f.write(resp.text)
            cif_path = f.name

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("X", cif_path)
        pdb_path = cif_path.replace(".cif", ".pdb")
        io_obj = PDBIO()
        io_obj.set_structure(structure)
        io_obj.save(pdb_path)
        os.unlink(cif_path)
        return pdb_path

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as f:
        f.write(resp.text)
        return f.name


def _run_propka_on_file(pdb_path: str, chain_id: str = "A") -> List[dict]:
    """
    Run PROPKA on a PDB file and return per-residue pKa results.
    Returns list of {resnum, resname, chain, pka}.
    """
    if not HAS_PROPKA:
        return _fallback_standard_pka(pdb_path, chain_id)

    try:
        # Run PROPKA
        mol = propka.run.single(pdb_path, optargs=["--quiet"])

        results = []
        for group in mol.conformations["AVR"].groups:
            resname = group.residue_type.strip().upper()
            if resname not in TITRATABLE_SET:
                continue
            if chain_id and group.atom.chain_id != chain_id:
                continue

            results.append({
                "resnum": group.atom.res_num,
                "resname": resname,
                "chain": group.atom.chain_id,
                "pka": round(group.pka_value, 2),
            })

        if results:
            return results

    except Exception as e:
        LOG.warning(f"PROPKA failed: {e}, falling back to standard pKa")

    return _fallback_standard_pka(pdb_path, chain_id)


def _fallback_standard_pka(pdb_path: str, chain_id: str) -> List[dict]:
    """Parse structure and assign standard textbook pKa values."""
    results = []
    if not HAS_BIOPYTHON:
        return results

    try:
        if pdb_path.lower().endswith(".cif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure("X", pdb_path)
        for model in structure:
            for chain in model:
                if chain.id != chain_id:
                    continue
                for residue in chain:
                    if not is_aa(residue, standard=True):
                        continue
                    resname = residue.get_resname().upper()
                    if resname in TITRATABLE_SET:
                        results.append({
                            "resnum": residue.id[1],
                            "resname": resname,
                            "chain": chain.id,
                            "pka": STANDARD_PKA[resname],
                        })
            break  # first model only
    except Exception as e:
        LOG.warning(f"Fallback parsing failed: {e}")

    return results


def _compute_pI(titratable_residues: List[dict]) -> float:
    """Binary search for the isoelectric point."""
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        charge = 0.0
        for r in titratable_residues:
            pka = r["pka"]
            rtype = TITRATABLE_TYPE.get(r["resname"], "acid")
            if rtype == "acid":
                charge += -1.0 / (1.0 + 10 ** (pka - mid))
            else:
                charge += 1.0 / (1.0 + 10 ** (mid - pka))
        if charge > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


# ═══════════════════════════════════════════════════════════════
# API Endpoints
# ═══════════════════════════════════════════════════════════════

@app.on_event("startup")
async def startup():
    """Load BRENDA query cache and credentials on server start."""
    _load_query_cache()
    _load_brenda_credentials()


@app.get("/health", response_model=HealthResponse)
async def health_check():
    """Check server status and available features."""
    creds = _load_brenda_credentials()
    return HealthResponse(
        status="ok",
        propka_available=HAS_PROPKA,
        biopython_available=HAS_BIOPYTHON,
        brenda_soap_available=HAS_ZEEP,
        brenda_credentials_set=creds is not None,
        brenda_cached_ec_count=len(_brenda_query_cache),
    )


@app.get("/propka/{pdb_id}", response_model=PropkaResult)
async def analyze_pdb(pdb_id: str, chain: str = "A"):
    """
    Run PROPKA analysis on a PDB structure from RCSB.

    Example: GET /propka/1LYZ?chain=A
    """
    pdb_path = None
    try:
        pdb_path = _download_pdb_file(pdb_id)
        raw_results = _run_propka_on_file(pdb_path, chain)

        titratable = []
        for r in raw_results:
            resname = r["resname"]
            std = STANDARD_PKA.get(resname, r["pka"])
            titratable.append(PropkaResidue(
                pos=r["resnum"],
                one=THREE_TO_ONE.get(resname, "X"),
                three=resname,
                chain=r["chain"],
                pka=r["pka"],
                type=TITRATABLE_TYPE.get(resname, "acid"),
                cls=CLS_MAP.get(resname, "Other"),
                std_pka=std,
                pka_shift=round(r["pka"] - std, 2),
            ))

        pI = _compute_pI(raw_results) if raw_results else 7.0

        return PropkaResult(
            pdb_id=pdb_id.upper(),
            chain_id=chain,
            n_residues=0,
            n_titratable=len(titratable),
            titratable=titratable,
            pI=pI,
            method="propka" if HAS_PROPKA else "standard_fallback",
        )

    finally:
        if pdb_path and os.path.exists(pdb_path):
            try:
                os.unlink(pdb_path)
            except Exception:
                pass


@app.post("/propka/upload", response_model=PropkaResult)
async def analyze_upload(
    file: UploadFile = File(...),
    chain: str = "A",
):
    """
    Run PROPKA analysis on an uploaded PDB/mmCIF file.

    Example: POST /propka/upload  (multipart form with 'file' field)
    """
    content = await file.read()
    filename = file.filename or "upload.pdb"

    suffix = ".cif" if filename.lower().endswith((".cif", ".mmcif")) else ".pdb"

    pdb_path = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=suffix, delete=False, mode="wb"
        ) as f:
            f.write(content)
            raw_path = f.name

        # If mmCIF, convert to PDB for PROPKA
        if suffix == ".cif" and HAS_BIOPYTHON:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("X", raw_path)
            pdb_path = raw_path.replace(".cif", ".pdb")
            io_obj = PDBIO()
            io_obj.set_structure(structure)
            io_obj.save(pdb_path)
            os.unlink(raw_path)
        else:
            pdb_path = raw_path

        raw_results = _run_propka_on_file(pdb_path, chain)

        titratable = []
        for r in raw_results:
            resname = r["resname"]
            std = STANDARD_PKA.get(resname, r["pka"])
            titratable.append(PropkaResidue(
                pos=r["resnum"],
                one=THREE_TO_ONE.get(resname, "X"),
                three=resname,
                chain=r["chain"],
                pka=r["pka"],
                type=TITRATABLE_TYPE.get(resname, "acid"),
                cls=CLS_MAP.get(resname, "Other"),
                std_pka=std,
                pka_shift=round(r["pka"] - std, 2),
            ))

        pI = _compute_pI(raw_results) if raw_results else 7.0

        # Extract PDB ID from file if possible
        pdb_id = "UPLOAD"
        text = content.decode("utf-8", errors="replace")
        for line in text.split("\n")[:100]:
            if line.startswith("HEADER") and len(line) >= 66:
                pdb_id = line[62:66].strip() or "UPLOAD"
                break
            if line.startswith("_entry.id"):
                parts = line.split()
                if len(parts) >= 2:
                    pdb_id = parts[1].strip().upper()
                    break

        return PropkaResult(
            pdb_id=pdb_id,
            chain_id=chain,
            n_residues=0,
            n_titratable=len(titratable),
            titratable=titratable,
            pI=pI,
            method="propka" if HAS_PROPKA else "standard_fallback",
        )

    finally:
        if pdb_path and os.path.exists(pdb_path):
            try:
                os.unlink(pdb_path)
            except Exception:
                pass


@app.get("/brenda/{ec_number}", response_model=Optional[BrendaResult])
async def lookup_brenda_ph(ec_number: str):
    """
    Look up pH optimum data from BRENDA SOAP API for a given EC number.

    Requires BRENDA credentials to be configured first via POST /brenda/credentials.
    Results are cached locally to minimize API calls.

    Example: GET /brenda/3.2.1.17
    """
    entries = _query_brenda_ph_optimum(ec_number)

    if not entries:
        raise HTTPException(404, f"No BRENDA pH data found for EC {ec_number}")

    ph_values = [e["ph"] for e in entries]
    organisms = list(set(e.get("organism", "unknown") for e in entries))

    return BrendaResult(
        ec_number=ec_number,
        ph_values=ph_values,
        ph_mean=round(float(np.mean(ph_values)), 2),
        ph_median=round(float(np.median(ph_values)), 2),
        ph_min=round(min(ph_values), 2),
        ph_max=round(max(ph_values), 2),
        n_entries=len(entries),
        organisms=organisms[:20],
    )


# —— Pydantic models for credential management ——

class BrendaCredentialsInput(BaseModel):
    email: str
    password: str  # plain text — hashed before storage


class BrendaCredentialsStatus(BaseModel):
    configured: bool
    email: Optional[str] = None
    message: str


@app.post("/brenda/credentials", response_model=BrendaCredentialsStatus)
async def set_brenda_credentials(creds: BrendaCredentialsInput):
    """
    Set BRENDA API credentials (Bring Your Own Key).

    Register for a free account at https://www.brenda-enzymes.org/register.php
    then provide your email and password here. The password is stored locally
    as a SHA-256 hash (the format BRENDA's SOAP API requires).

    This only needs to be done once — credentials persist across sessions.
    """
    if not creds.email or not creds.password:
        raise HTTPException(400, "Both email and password are required")

    _save_brenda_credentials(creds.email, creds.password)

    return BrendaCredentialsStatus(
        configured=True,
        email=creds.email,
        message="BRENDA credentials saved successfully. You can now query pH data.",
    )


@app.get("/brenda/credentials/status", response_model=BrendaCredentialsStatus)
async def check_brenda_credentials():
    """Check if BRENDA credentials are configured."""
    creds = _load_brenda_credentials()
    if creds:
        return BrendaCredentialsStatus(
            configured=True,
            email=creds["email"],
            message="BRENDA credentials are configured.",
        )
    return BrendaCredentialsStatus(
        configured=False,
        message="BRENDA credentials not set. Use POST /brenda/credentials "
                "with your email and password. Register free at "
                "https://www.brenda-enzymes.org/register.php",
    )


@app.delete("/brenda/credentials")
async def delete_brenda_credentials():
    """Remove stored BRENDA credentials."""
    global _brenda_credentials, _brenda_soap_client
    _brenda_credentials = None
    _brenda_soap_client = None
    if BRENDA_CREDS_FILE.exists():
        BRENDA_CREDS_FILE.unlink()
    return {"status": "deleted", "message": "BRENDA credentials removed."}


@app.get("/brenda/stats")
async def brenda_stats():
    """Return stats about the cached BRENDA queries."""
    creds = _load_brenda_credentials()
    return {
        "soap_available": HAS_ZEEP,
        "credentials_configured": creds is not None,
        "cached_ec_numbers": len(_brenda_query_cache),
        "total_cached_ph_entries": sum(len(v) for v in _brenda_query_cache.values()),
    }


@app.post("/brenda/cache/clear")
async def brenda_clear_cache():
    """Clear the local BRENDA query cache (forces fresh API calls)."""
    global _brenda_query_cache
    _brenda_query_cache = {}
    if BRENDA_QUERY_CACHE_FILE.exists():
        BRENDA_QUERY_CACHE_FILE.unlink()
        LOG.info("Deleted BRENDA query cache")
    return {
        "status": "cleared",
        "message": "BRENDA query cache cleared. Next queries will hit the SOAP API.",
    }


# ═══════════════════════════════════════════════════════════════
# Run with: python conformseek_server.py
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import uvicorn
    creds = _load_brenda_credentials()
    print("=" * 60)
    print("  ConformSeek API Server v2.0 (BRENDA SOAP API)")
    print("=" * 60)
    print(f"  PROPKA:    {'available' if HAS_PROPKA else 'NOT INSTALLED (pip install propka)'}")
    print(f"  BioPython: {'available' if HAS_BIOPYTHON else 'NOT INSTALLED (pip install biopython)'}")
    print(f"  Zeep SOAP: {'available' if HAS_ZEEP else 'NOT INSTALLED (pip install zeep)'}")
    print(f"  BRENDA:    {'credentials set (' + creds['email'] + ')' if creds else 'NOT CONFIGURED'}")
    if not creds:
        print(f"             POST /brenda/credentials to set your BRENDA login")
        print(f"             Register free: https://www.brenda-enzymes.org/register.php")
    print(f"  Cache:     {len(_brenda_query_cache)} EC numbers cached")
    print(f"  Docs:      http://localhost:8000/docs")
    print("=" * 60)
    uvicorn.run(app, host="0.0.0.0", port=8000)
