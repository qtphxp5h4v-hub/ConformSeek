/* global PDBeMolstarPlugin */
import { useState, useEffect, useRef, useCallback, forwardRef, useImperativeHandle } from "react";
import * as d3 from "d3";

// ═══════════════════════════════════════════════════════════════
// Constants
// ═══════════════════════════════════════════════════════════════
const THREE_TO_FULL = {
  ALA:"Alanine",ARG:"Arginine",ASN:"Asparagine",ASP:"Aspartic Acid",
  CYS:"Cysteine",GLN:"Glutamine",GLU:"Glutamic Acid",GLY:"Glycine",
  HIS:"Histidine",ILE:"Isoleucine",LEU:"Leucine",LYS:"Lysine",
  MET:"Methionine",PHE:"Phenylalanine",PRO:"Proline",SER:"Serine",
  THR:"Threonine",TRP:"Tryptophan",TYR:"Tyrosine",VAL:"Valine"
};
const ONE_TO_THREE = {A:"ALA",R:"ARG",N:"ASN",D:"ASP",C:"CYS",Q:"GLN",E:"GLU",G:"GLY",H:"HIS",I:"ILE",L:"LEU",K:"LYS",M:"MET",F:"PHE",P:"PRO",S:"SER",T:"THR",W:"TRP",Y:"TYR",V:"VAL"};
const TITRATABLE = {ASP:{cls:"Acidic",stdPka:3.65},GLU:{cls:"Acidic",stdPka:4.25},HIS:{cls:"Basic",stdPka:6.00},CYS:{cls:"Thiol",stdPka:8.18},LYS:{cls:"Basic",stdPka:10.53},TYR:{cls:"Phenolic",stdPka:10.07}};
const TITRATABLE_SET = new Set(Object.keys(TITRATABLE));
const CLS_COLORS = {Acidic:"#ff6b6b",Basic:"#4dabf7",Thiol:"#ffd43b",Phenolic:"#da77f2"};

// ═══════════════════════════════════════════════════════════════
// Backend API (ConformSeek server)
// ═══════════════════════════════════════════════════════════════
const API_BASE = "http://localhost:8000";

async function fetchPropkaFromBackend(pdbId, chain = "A") {
  try {
    const res = await fetch(`${API_BASE}/propka/${pdbId}?chain=${chain}`, { signal: AbortSignal.timeout(30000) });
    if (res.ok) {
      const data = await res.json();
      if (data.titratable && data.titratable.length > 0) {
        return data;
      }
    }
  } catch (e) {
    console.warn("PROPKA backend unavailable:", e.message || e);
  }
  return null;
}

async function checkBackendHealth() {
  try {
    const res = await fetch(`${API_BASE}/health`, { signal: AbortSignal.timeout(3000) });
    if (res.ok) return await res.json();
  } catch (e) { /* backend not running */ }
  return null;
}

async function fetchBrendaPH(ecNumber) {
  if (!ecNumber) return null;
  try {
    const res = await fetch(`${API_BASE}/brenda/${ecNumber}`, { signal: AbortSignal.timeout(10000) });
    if (res.ok) return await res.json(); // { ec_number, ph_mean, ph_median, ph_min, ph_max, n_entries, organisms }
    if (res.status === 401) return { _needsCredentials: true };
  } catch (e) {
    console.warn("BRENDA lookup failed:", e.message || e);
  }
  return null;
}

// ── BRENDA Credential Management ──
async function fetchBrendaCredentialStatus() {
  try {
    const res = await fetch(`${API_BASE}/brenda/credentials/status`, { signal: AbortSignal.timeout(3000) });
    if (res.ok) return await res.json();
  } catch (e) { /* backend not running */ }
  return null;
}

async function saveBrendaCredentials(email, password) {
  try {
    const res = await fetch(`${API_BASE}/brenda/credentials`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ email, password }),
      signal: AbortSignal.timeout(5000),
    });
    return await res.json();
  } catch (e) {
    console.warn("Failed to save BRENDA credentials:", e.message || e);
    return null;
  }
}

async function deleteBrendaCredentials() {
  try {
    const res = await fetch(`${API_BASE}/brenda/credentials`, {
      method: "DELETE",
      signal: AbortSignal.timeout(3000),
    });
    return await res.json();
  } catch (e) {
    console.warn("Failed to delete BRENDA credentials:", e.message || e);
    return null;
  }
}

async function clearBrendaCache() {
  try {
    const res = await fetch(`${API_BASE}/brenda/cache/clear`, {
      method: "POST",
      signal: AbortSignal.timeout(3000),
    });
    return await res.json();
  } catch (e) {
    return null;
  }
}

// ═══════════════════════════════════════════════════════════════
// Subcellular localization → functional pH mapping (Tier 2)
// ═══════════════════════════════════════════════════════════════
// Priority tiers for subcellular matching (higher = preferred as primary function site)
// Tier 3: Primary functional compartments (where the protein does its job)
// Tier 2: Common compartments (reasonable default)
// Tier 1: Transient/recycling/processing compartments (less likely primary function)
const SUBCELLULAR_PH_MAP = {
  // Tier 3 — primary functional compartments
  "hemoglobin complex":    { ph: 7.4, range: [7.35, 7.45], label: "Blood (hemoglobin)", priority: 30 },
  "haptoglobin-hemoglobin complex": { ph: 7.4, range: [7.35, 7.45], label: "Blood (hemoglobin)", priority: 30 },
  "blood microparticle":   { ph: 7.4, range: [7.35, 7.45], label: "Blood", priority: 28 },
  "blood":                 { ph: 7.4, range: [7.35, 7.45], label: "Blood", priority: 28 },
  "secreted":              { ph: 7.4, range: [7.35, 7.45], label: "Blood/secreted", priority: 25 },
  "lysosome":              { ph: 4.7, range: [4.5, 5.0],  label: "Lysosome", priority: 25 },
  "lysosomal lumen":       { ph: 4.7, range: [4.5, 5.0],  label: "Lysosomal lumen", priority: 25 },
  "thylakoid lumen":       { ph: 5.7, range: [5.0, 6.5],  label: "Thylakoid lumen", priority: 25 },
  "chloroplast stroma":    { ph: 8.0, range: [7.5, 8.5],  label: "Chloroplast stroma", priority: 25 },
  "mitochondrial matrix":  { ph: 7.8, range: [7.5, 8.0],  label: "Mitochondrial matrix", priority: 24 },
  "mitochondrial intermembrane space": { ph: 6.8, range: [6.5, 7.0], label: "Mito IMS", priority: 24 },
  "stomach":               { ph: 2.0, range: [1.5, 3.5],  label: "Gastric lumen", priority: 25 },
  "vacuole":               { ph: 5.5, range: [5.0, 6.0],  label: "Vacuole", priority: 24 },
  "apoplast":              { ph: 5.5, range: [5.0, 6.5],  label: "Apoplast", priority: 24 },
  "peroxisome":            { ph: 7.0, range: [6.9, 7.1],  label: "Peroxisome", priority: 23 },

  // Tier 2 — common compartments
  "nucleus":               { ph: 7.2, range: [7.0, 7.4],  label: "Nucleus", priority: 15 },
  "nucleoplasm":           { ph: 7.2, range: [7.0, 7.4],  label: "Nucleoplasm", priority: 15 },
  "cytoplasm":             { ph: 7.2, range: [7.0, 7.4],  label: "Cytoplasm", priority: 14 },
  "cytosol":               { ph: 7.2, range: [7.0, 7.4],  label: "Cytosol", priority: 14 },
  "mitochondrion":         { ph: 7.8, range: [7.5, 8.0],  label: "Mitochondrial matrix", priority: 14 },
  "extracellular region":  { ph: 7.4, range: [7.2, 7.5],  label: "Extracellular", priority: 13 },
  "extracellular space":   { ph: 7.4, range: [7.2, 7.5],  label: "Extracellular", priority: 13 },
  "extracellular":         { ph: 7.4, range: [7.2, 7.5],  label: "Extracellular", priority: 13 },
  "cell surface":          { ph: 7.4, range: [7.2, 7.5],  label: "Cell surface", priority: 13 },
  "extracellular matrix":  { ph: 7.4, range: [7.2, 7.5],  label: "Extracellular matrix", priority: 13 },
  "plasma membrane":       { ph: 7.2, range: [7.0, 7.4],  label: "Plasma membrane", priority: 12 },
  "endoplasmic reticulum": { ph: 7.2, range: [7.0, 7.4],  label: "ER", priority: 12 },
  "endoplasmic reticulum lumen": { ph: 7.2, range: [7.0, 7.4], label: "ER lumen", priority: 12 },
  "golgi apparatus":       { ph: 6.4, range: [6.0, 6.7],  label: "Golgi apparatus", priority: 12 },
  "cis-golgi":             { ph: 6.7, range: [6.5, 7.0],  label: "cis-Golgi", priority: 12 },
  "trans-golgi":           { ph: 6.2, range: [6.0, 6.4],  label: "trans-Golgi", priority: 12 },
  "intestine":             { ph: 7.0, range: [6.0, 8.0],  label: "Intestinal lumen", priority: 12 },

  // Tier 1 — transient/recycling/processing (deprioritized)
  "endosome":              { ph: 5.5, range: [5.0, 6.5],  label: "Endosome", priority: 5 },
  "early endosome":        { ph: 6.2, range: [5.9, 6.5],  label: "Early endosome", priority: 5 },
  "late endosome":         { ph: 5.5, range: [5.0, 6.0],  label: "Late endosome", priority: 5 },
  "recycling endosome":    { ph: 6.4, range: [6.2, 6.5],  label: "Recycling endosome", priority: 5 },
  "endocytic vesicle lumen": { ph: 5.5, range: [5.0, 6.5], label: "Endocytic vesicle", priority: 3 },
  "extracellular exosome": { ph: 7.4, range: [7.2, 7.5],  label: "Extracellular exosome", priority: 3 },
  "membrane":              { ph: 7.2, range: [7.0, 7.4],  label: "Membrane", priority: 2 },
};

async function fetchSubcellularPH(uniprotAcc) {
  if (!uniprotAcc) return null;
  try {
    const res = await fetch(`https://rest.uniprot.org/uniprotkb/${uniprotAcc}.json`, {
      signal: AbortSignal.timeout(8000),
      headers: { Accept: "application/json" },
    });
    if (!res.ok) return null;
    const data = await res.json();

    // Collect location terms from multiple sources
    const locations = [];

    // Source 1: SUBCELLULAR LOCATION comments
    for (const c of (data?.comments || [])) {
      if (c.commentType === "SUBCELLULAR LOCATION") {
        for (const sl of (c.subcellularLocations || [])) {
          if (sl.location?.value) locations.push(sl.location.value.toLowerCase());
          if (sl.topology?.value) locations.push(sl.topology.value.toLowerCase());
        }
      }
    }

    // Source 2: Keywords (Secreted, Lysosome, Nucleus, etc.)
    for (const kw of (data?.keywords || [])) {
      const v = (kw.name || "").toLowerCase();
      if (["secreted", "lysosome", "nucleus", "mitochondrion", "membrane",
           "cytoplasm", "endoplasmic reticulum", "golgi apparatus", "peroxisome",
           "extracellular matrix", "cell membrane"].includes(v)) {
        locations.push(v);
      }
    }

    // Source 3: GO cellular component terms (cross-references)
    for (const xref of (data?.uniProtKBCrossReferences || [])) {
      if (xref.database === "GO") {
        const props = xref.properties || [];
        const term = props.find(p => p.key === "GoTerm")?.value || "";
        if (term.startsWith("C:")) {
          locations.push(term.substring(2).toLowerCase());
        }
      }
    }

    // Source 4: FUNCTION comment may mention blood, plasma, etc.
    for (const c of (data?.comments || [])) {
      if (c.commentType === "FUNCTION") {
        const txt = (c.texts?.[0]?.value || "").toLowerCase();
        if (txt.includes("blood") || txt.includes("plasma") || txt.includes("serum")) locations.push("blood");
        if (txt.includes("oxygen transport")) locations.push("blood");
        if (txt.includes("gastric") || txt.includes("stomach")) locations.push("stomach");
        if (txt.includes("lysosom")) locations.push("lysosome");
      }
    }

    console.log("UniProt locations found:", locations);

    if (locations.length === 0) return null;

    // Match — prefer highest priority (primary functional compartment)
    let bestMatch = null;
    let bestPriority = -1;
    for (const loc of locations) {
      for (const [key, mapping] of Object.entries(SUBCELLULAR_PH_MAP)) {
        if (loc.includes(key) && mapping.priority > bestPriority) {
          bestPriority = mapping.priority;
          bestMatch = { ...mapping, rawLocation: loc, source: "UniProt subcellular localization" };
        }
      }
    }
    return bestMatch ? { ...bestMatch, allLocations: [...new Set(locations)] } : null;
  } catch (e) {
    console.warn("UniProt subcellular lookup failed:", e.message || e);
    return null;
  }
}

function propkaToTitratable(propkaData) {
  return propkaData.titratable.map(r => ({
    pos: r.pos, one: r.one, three: r.three,
    full: THREE_TO_FULL[r.three] || r.three,
    cls: r.cls, pka: r.pka,
    type: r.type,
    stdPka: r.std_pka,
    pkaShift: r.pka_shift,
  }));
}

// ═══════════════════════════════════════════════════════════════
// Titration computation (Henderson-Hasselbalch)
// ═══════════════════════════════════════════════════════════════
function computeCharge(titratable, ph) {
  let charge = 0;
  for (const r of titratable) {
    const pka = r.pka;
    if (r.type === "acid") {
      charge += -1 / (1 + Math.pow(10, pka - ph));
    } else {
      charge += 1 / (1 + Math.pow(10, ph - pka));
    }
  }
  return charge;
}

function computeTitrationCurve(titratable) {
  const points = [];
  for (let ph = 0; ph <= 14; ph += 0.1) {
    points.push({ ph: Math.round(ph * 10) / 10, charge: computeCharge(titratable, ph) });
  }
  return points;
}

function computePI(titratable) {
  let lo = 0, hi = 14;
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2;
    const c = computeCharge(titratable, mid);
    if (c > 0) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

// ═══════════════════════════════════════════════════════════════
// UniProt → PDB ID mapping via RCSB Search API
// ═══════════════════════════════════════════════════════════════
async function lookupUniProt(uniprotId) {
  // First try to find PDB structures mapped to this UniProt ID
  try {
    const query = {
      query: {
        type: "terminal",
        service: "text",
        parameters: {
          attribute: "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
          operator: "exact_match",
          value: uniprotId.toUpperCase(),
        },
      },
      request_options: { paginate: { start: 0, rows: 10 } },
      return_type: "entry",
    };
    const res = await fetch(`https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify(query))}`);
    if (res.ok) {
      const data = await res.json();
      const pdbIds = data?.result_set?.map(r => r.identifier) || [];
      if (pdbIds.length > 0) {
        return { type: "pdb_list", pdbIds, uniprotId };
      }
    }
  } catch (e) {
    console.error("RCSB search failed:", e);
  }

  // If no PDB structures, try AlphaFold
  try {
    const afRes = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`);
    if (afRes.ok) {
      const afData = await afRes.json();
      if (afData.length > 0) {
        return {
          type: "alphafold",
          url: afData[0].cifUrl,
          entryId: afData[0].entryId,
          uniprotId,
        };
      }
    }
  } catch (e) {
    console.error("AlphaFold lookup failed:", e);
  }

  return { type: "not_found", uniprotId };
}

// ═══════════════════════════════════════════════════════════════
// Parse sequence from uploaded PDB/CIF file text
// ═══════════════════════════════════════════════════════════════
function parseSequenceFromFileText(text, format) {
  const THREE_TO_ONE_MAP = {
    ALA:"A",ARG:"R",ASN:"N",ASP:"D",CYS:"C",GLN:"Q",GLU:"E",GLY:"G",
    HIS:"H",ILE:"I",LEU:"L",LYS:"K",MET:"M",PHE:"F",PRO:"P",SER:"S",
    THR:"T",TRP:"W",TYR:"Y",VAL:"V"
  };

  const residues = new Map(); // resNum -> resName

  if (format === "pdb") {
    for (const line of text.split("\n")) {
      if (line.startsWith("ATOM") || line.startsWith("HETATM")) {
        const chain = line[21];
        if (chain !== "A" && chain !== " ") continue; // default to chain A
        const resName = line.substring(17, 20).trim();
        const resNum = parseInt(line.substring(22, 26).trim());
        if (THREE_TO_ONE_MAP[resName] && !residues.has(resNum)) {
          residues.set(resNum, resName);
        }
      }
    }
  } else {
    // mmCIF - extract from _atom_site records
    for (const line of text.split("\n")) {
      if (line.startsWith("ATOM") || line.startsWith("HETATM")) {
        const parts = line.split(/\s+/);
        // mmCIF ATOM columns: group_PDB id type_symbol label_atom_id ... label_comp_id label_asym_id ... label_seq_id
        if (parts.length >= 9) {
          const resName = parts[5];
          const resNum = parseInt(parts[8]);
          if (THREE_TO_ONE_MAP[resName] && !isNaN(resNum) && !residues.has(resNum)) {
            residues.set(resNum, resName);
          }
        }
      }
    }
  }

  const sorted = [...residues.entries()].sort((a, b) => a[0] - b[0]);
  return sorted.map(([, name]) => THREE_TO_ONE_MAP[name] || "X").join("");
}

// ═══════════════════════════════════════════════════════════════
// Parse mmCIF/PDB sequence + metadata from RCSB API
// ═══════════════════════════════════════════════════════════════
async function fetchPDBInfo(pdbId) {
  const id = pdbId.toLowerCase();
  try {
    const res = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${id}`);
    const data = await res.json();
    const title = data?.struct?.title || "";
    const resolution = data?.rcsb_entry_info?.resolution_combined?.[0] || null;
    const method = data?.exptl?.[0]?.method || "";
    // Crystal growth conditions
    const crystGrow = data?.exptl_crystal_grow?.[0] || {};
    const crystPh = crystGrow?.pH || null;
    const crystTemp = crystGrow?.temp || null;
    const crystMethod = crystGrow?.method || null;
    const crystDetails = crystGrow?.pdbx_details || null;
    // Data collection temperature
    const diffrn = data?.diffrn?.[0] || {};
    const collectionTemp = diffrn?.ambient_temp || null;
    // Refinement stats
    const refine = data?.refine?.[0] || {};
    const rFree = refine?.ls_R_factor_R_free || null;
    const rWork = refine?.ls_R_factor_R_work || null;
    // Deposition info
    const deposited = data?.rcsb_accession_info?.deposit_date || null;
    const released = data?.rcsb_accession_info?.initial_release_date || null;
    // Organism
    const organism = data?.rcsb_entry_info?.polymer_entity_count_protein
      ? null : null; // will fetch from entity
    // Keywords
    const keywords = data?.struct_keywords?.pdbx_keywords || null;

    // Get polymer entity info for sequence + organism
    const polyRes = await fetch(`https://data.rcsb.org/rest/v1/core/polymer_entity/${id}/1`);
    const polyData = await polyRes.json();
    const sequence = polyData?.entity_poly?.pdbx_seq_one_letter_code_can || "";
    const orgName = polyData?.rcsb_entity_source_organism?.[0]?.ncbi_scientific_name || null;
    const expressionSystem = polyData?.rcsb_entity_source_organism?.[0]?.expression_system_ncbi_scientific_name || null;
    // EC number from enzyme classification (try multiple paths)
    let ecNumber = null;
    // Path 1: polymer entity rcsb_ec_lineage (full 4-level EC numbers)
    const ecLineage = polyData?.rcsb_polymer_entity?.rcsb_ec_lineage || [];
    const fullEcs = ecLineage.filter(e => e?.id && e.id.split(".").length === 4).map(e => e.id);
    if (fullEcs.length > 0) ecNumber = fullEcs[0];
    // Path 2: enzyme classification at entity level
    if (!ecNumber) {
      const enzClass = polyData?.rcsb_polymer_entity?.rcsb_enzyme_class_combined || [];
      const fullEcs2 = enzClass.filter(e => e?.ec && e.ec.split(".").length === 4).map(e => e.ec);
      if (fullEcs2.length > 0) ecNumber = fullEcs2[0];
    }
    // Path 3: pdbx_description might contain EC info
    if (!ecNumber) {
      const desc = polyData?.rcsb_polymer_entity?.pdbx_description || "";
      const ecMatch = desc.match(/EC\s*([\d]+\.[\d]+\.[\d]+\.[\d]+)/i);
      if (ecMatch) ecNumber = ecMatch[1];
    }
    // UniProt accession for subcellular localization lookup
    const refSeqIds = polyData?.rcsb_polymer_entity_container_identifiers?.reference_sequence_identifiers || [];
    const uniprotAcc = refSeqIds.find(r => r.database_name === "UniProt")?.database_accession || null;

    return {
      title, resolution, method, crystPh, crystTemp, crystMethod, crystDetails,
      collectionTemp, rFree, rWork, deposited, released, keywords,
      sequence, organism: orgName, expressionSystem, ecNumber, uniprotAcc
    };
  } catch (e) {
    console.error("Failed to fetch PDB info:", e);
    return {
      title: "", resolution: null, method: "", crystPh: null, crystTemp: null,
      crystMethod: null, crystDetails: null, collectionTemp: null, rFree: null,
      rWork: null, deposited: null, released: null, keywords: null,
      sequence: "", organism: null, expressionSystem: null, ecNumber: null, uniprotAcc: null
    };
  }
}

// ═══════════════════════════════════════════════════════════════
// Build titratable residue list from sequence
// ═══════════════════════════════════════════════════════════════
function buildTitratableFromSequence(sequence) {
  const residues = [];
  for (let i = 0; i < sequence.length; i++) {
    const one = sequence[i];
    const three = ONE_TO_THREE[one];
    if (three && TITRATABLE_SET.has(three)) {
      const info = TITRATABLE[three];
      // Add slight random perturbation to mimic environment effects
      // In real app this would come from PROPKA
      const envShift = (Math.random() - 0.5) * 1.5;
      residues.push({
        pos: i + 1,
        one,
        three,
        full: THREE_TO_FULL[three],
        cls: info.cls,
        pka: Math.round((info.stdPka + envShift) * 100) / 100,
        type: info.cls === "Basic" ? "base" : "acid",
      });
    }
  }
  return residues;
}

// ═══════════════════════════════════════════════════════════════
// TitrationCurve component with embedded pH line
// ═══════════════════════════════════════════════════════════════
function TitrationCurve({ curve, ph, pI, crystPh, brendaPh, onPhChange }) {
  const svgRef = useRef();
  const margin = { top: 20, right: 30, bottom: 40, left: 55 };
  const w = 480, h = 240;
  const innerW = w - margin.left - margin.right;
  const innerH = h - margin.top - margin.bottom;

  useEffect(() => {
    if (!curve || curve.length === 0) return;
    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    const x = d3.scaleLinear().domain([0, 14]).range([0, innerW]);
    const yExt = d3.extent(curve, d => d.charge);
    const yPad = (yExt[1] - yExt[0]) * 0.1;
    const y = d3.scaleLinear().domain([yExt[0] - yPad, yExt[1] + yPad]).range([innerH, 0]);

    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    // Grid
    g.append("g").attr("transform", `translate(0,${innerH})`).call(d3.axisBottom(x).ticks(14).tickSize(-innerH)).selectAll("line").attr("stroke", "#1a2030").attr("stroke-opacity", 0.6);
    g.append("g").call(d3.axisLeft(y).ticks(6).tickSize(-innerW)).selectAll("line").attr("stroke", "#1a2030").attr("stroke-opacity", 0.6);
    g.selectAll(".domain").attr("stroke", "#334");
    g.selectAll("text").attr("fill", "#8899aa").attr("font-size", "10px");

    // Zero line
    g.append("line").attr("x1", 0).attr("x2", innerW).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", "#445").attr("stroke-dasharray", "4,3");

    // pH-sensitive band
    const bandLeft = Math.max(0, ph - 1);
    const bandRight = Math.min(14, ph + 1);
    g.append("rect").attr("x", x(bandLeft)).attr("width", x(bandRight) - x(bandLeft)).attr("y", 0).attr("height", innerH).attr("fill", "#ff4466").attr("opacity", 0.06);

    // Titration curve
    const line = d3.line().x(d => x(d.ph)).y(d => y(d.charge)).curve(d3.curveCatmullRom);
    g.append("path").datum(curve).attr("d", line).attr("fill", "none").attr("stroke", "#4488ff").attr("stroke-width", 2.5);

    // pI line
    if (pI != null) {
      g.append("line").attr("x1", x(pI)).attr("x2", x(pI)).attr("y1", 0).attr("y2", innerH).attr("stroke", "#00cc88").attr("stroke-width", 1.5).attr("stroke-dasharray", "5,3");
      g.append("text").attr("x", x(pI) + 4).attr("y", 14).text(`pI ${pI.toFixed(1)}`).attr("fill", "#00cc88").attr("font-size", "10px");
    }

    // Cryst pH
    if (crystPh != null) {
      g.append("line").attr("x1", x(crystPh)).attr("x2", x(crystPh)).attr("y1", 0).attr("y2", innerH).attr("stroke", "#ffcc00").attr("stroke-width", 1.2).attr("stroke-dasharray", "3,4");
      g.append("text").attr("x", x(crystPh) + 4).attr("y", innerH - 6).text(`Cryst ${crystPh}`).attr("fill", "#ffcc00").attr("font-size", "9px");
    }

    // BRENDA functional pH (green dashed)
    if (brendaPh != null) {
      g.append("line").attr("x1", x(brendaPh)).attr("x2", x(brendaPh)).attr("y1", 0).attr("y2", innerH).attr("stroke", "#00ff88").attr("stroke-width", 1.5).attr("stroke-dasharray", "6,2");
      g.append("text").attr("x", x(brendaPh) + 4).attr("y", 28).text(`Func ${brendaPh.toFixed(1)}`).attr("fill", "#00ff88").attr("font-size", "9px").attr("font-weight", "600");
    }

    // Current pH indicator
    const currentCharge = computeCharge([], ph); // approximate from curve
    const closestPt = curve.reduce((prev, curr) => Math.abs(curr.ph - ph) < Math.abs(prev.ph - ph) ? curr : prev);
    g.append("line").attr("x1", x(ph)).attr("x2", x(ph)).attr("y1", 0).attr("y2", innerH).attr("stroke", "#ff4466").attr("stroke-width", 2).attr("stroke-dasharray", "6,3");
    g.append("circle").attr("cx", x(ph)).attr("cy", y(closestPt.charge)).attr("r", 6).attr("fill", "#ff4466").attr("stroke", "#fff").attr("stroke-width", 1.5);
    g.append("text").attr("x", x(ph) + 8).attr("y", y(closestPt.charge) + 4).text(`pH ${ph.toFixed(1)}: ${closestPt.charge >= 0 ? "+" : ""}${closestPt.charge.toFixed(1)}`).attr("fill", "#ff4466").attr("font-size", "11px").attr("font-weight", "600");

    // Click-to-set-pH
    g.append("rect").attr("width", innerW).attr("height", innerH).attr("fill", "transparent").attr("cursor", "crosshair")
      .on("click", (event) => {
        const [mx] = d3.pointer(event);
        const newPh = Math.round(x.invert(mx) * 10) / 10;
        if (newPh >= 0 && newPh <= 14) onPhChange(newPh);
      })
      .on("mousemove", function(event) {
        const [mx] = d3.pointer(event);
        const hoverPh = x.invert(mx);
        const hPt = curve.reduce((p, c) => Math.abs(c.ph - hoverPh) < Math.abs(p.ph - hoverPh) ? c : p);
        d3.select(this.parentNode).selectAll(".hover-line").remove();
        d3.select(this.parentNode).selectAll(".hover-dot").remove();
        d3.select(this.parentNode).append("line").attr("class", "hover-line").attr("x1", mx).attr("x2", mx).attr("y1", 0).attr("y2", innerH).attr("stroke", "#ffffff20").attr("stroke-width", 1);
        d3.select(this.parentNode).append("circle").attr("class", "hover-dot").attr("cx", mx).attr("cy", y(hPt.charge)).attr("r", 3).attr("fill", "#fff").attr("opacity", 0.5);
      })
      .on("mouseleave", function() {
        d3.select(this.parentNode).selectAll(".hover-line").remove();
        d3.select(this.parentNode).selectAll(".hover-dot").remove();
      });

    // Axis labels
    svg.append("text").attr("x", w / 2).attr("y", h - 2).attr("text-anchor", "middle").attr("fill", "#8899aa").attr("font-size", "11px").text("pH");
    svg.append("text").attr("transform", `rotate(-90)`).attr("x", -h / 2).attr("y", 14).attr("text-anchor", "middle").attr("fill", "#8899aa").attr("font-size", "11px").text("Net Charge");
  }, [curve, ph, pI, crystPh, brendaPh]);

  return <svg ref={svgRef} width={w} height={h} />;
}

// ═══════════════════════════════════════════════════════════════
// pKa Bar Chart
// ═══════════════════════════════════════════════════════════════
function PkaBarChart({ titratable, ph }) {
  const svgRef = useRef();
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };
  const w = 480, h = 240;

  useEffect(() => {
    if (!titratable || titratable.length === 0) return;
    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();
    const innerW = w - margin.left - margin.right;
    const innerH = h - margin.top - margin.bottom;
    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    const x = d3.scaleBand().domain(titratable.map(d => d.pos)).range([0, innerW]).padding(0.3);
    const y = d3.scaleLinear().domain([0, 14]).range([innerH, 0]);

    g.append("g").attr("transform", `translate(0,${innerH})`).call(d3.axisBottom(x).tickValues(titratable.filter((_, i) => i % Math.ceil(titratable.length / 20) === 0).map(d => d.pos))).selectAll("text").attr("fill", "#8899aa").attr("font-size", "8px");
    g.append("g").call(d3.axisLeft(y).ticks(7)).selectAll("text").attr("fill", "#8899aa").attr("font-size", "10px");
    g.selectAll(".domain").attr("stroke", "#334");

    // pH band
    g.append("rect").attr("x", 0).attr("width", innerW).attr("y", y(ph + 1)).attr("height", y(ph - 1) - y(ph + 1)).attr("fill", "#ff4466").attr("opacity", 0.06);
    g.append("line").attr("x1", 0).attr("x2", innerW).attr("y1", y(ph)).attr("y2", y(ph)).attr("stroke", "#ff4466").attr("stroke-width", 1.5).attr("stroke-dasharray", "5,3");

    // Bars
    g.selectAll(".bar").data(titratable).enter().append("rect")
      .attr("x", d => x(d.pos))
      .attr("y", d => y(d.pka))
      .attr("width", x.bandwidth())
      .attr("height", d => innerH - y(d.pka))
      .attr("fill", d => {
        const delta = Math.abs(d.pka - ph);
        return delta < 1 ? "#ff4466" : delta < 2 ? "#ffcc00" : "#4488ff";
      })
      .attr("opacity", 0.85)
      .attr("rx", 1);

    svg.append("text").attr("x", w / 2).attr("y", h - 2).attr("text-anchor", "middle").attr("fill", "#8899aa").attr("font-size", "11px").text("Residue Position");
    svg.append("text").attr("transform", "rotate(-90)").attr("x", -h / 2).attr("y", 14).attr("text-anchor", "middle").attr("fill", "#8899aa").attr("font-size", "11px").text("pKa");
  }, [titratable, ph]);

  return <svg ref={svgRef} width={w} height={h} />;
}

// ═══════════════════════════════════════════════════════════════
// Titratable Residues Table
// ═══════════════════════════════════════════════════════════════
function TitratableTable({ titratable, ph, pkaSource }) {
  const [sortKey, setSortKey] = useState("pos");
  const [sortAsc, setSortAsc] = useState(true);
  const hasPropka = pkaSource === "propka" || pkaSource === "standard_fallback";

  const sorted = [...titratable].sort((a, b) => {
    const va = a[sortKey], vb = b[sortKey];
    return sortAsc ? (va < vb ? -1 : 1) : (va > vb ? -1 : 1);
  });

  const handleSort = (key) => {
    if (sortKey === key) setSortAsc(!sortAsc);
    else { setSortKey(key); setSortAsc(true); }
  };

  const arrow = (key) => sortKey === key ? (sortAsc ? " ↑" : " ↓") : "";

  const columns = [["pos","Pos"],["one","AA"],["full","Amino Acid"],["cls","Class"],["pka","pKa"]];
  if (hasPropka) columns.push(["stdPka","Std pKa"],["pkaShift","Env Shift"]);
  columns.push(["delta","ΔpKa"]);

  return (
    <div style={{ maxHeight: 360, overflowY: "auto", borderRadius: 8, border: "1px solid #1e2a3a" }}>
      <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 12 }}>
        <thead>
          <tr style={{ position: "sticky", top: 0, background: "#0d1520", zIndex: 1 }}>
            {columns.map(([k, label]) => (
              <th key={k} onClick={() => handleSort(k)} style={{ padding: "8px 10px", textAlign: "left", cursor: "pointer", color: "#8899aa", borderBottom: "1px solid #1e2a3a", whiteSpace: "nowrap", userSelect: "none" }}>
                {label}{arrow(k)}
              </th>
            ))}
            <th style={{ padding: "8px 10px", textAlign: "left", color: "#8899aa", borderBottom: "1px solid #1e2a3a" }}>State</th>
          </tr>
        </thead>
        <tbody>
          {sorted.map((r, i) => {
            const delta = r.pka - ph;
            const protonated = r.type === "acid" ? (ph < r.pka) : (ph < r.pka);
            const state = protonated ? "Protonated" : "Deprotonated";
            const sensitive = Math.abs(delta) < 1;
            return (
              <tr key={r.pos} style={{ background: i % 2 === 0 ? "#0a1018" : "#0d1520", borderBottom: "1px solid #111a25" }}>
                <td style={{ padding: "6px 10px", color: "#ccd" }}>{r.pos}</td>
                <td style={{ padding: "6px 10px", color: CLS_COLORS[r.cls] || "#ccd", fontWeight: 700, fontFamily: "monospace" }}>{r.one}</td>
                <td style={{ padding: "6px 10px", color: "#ccd" }}>{r.full}</td>
                <td style={{ padding: "6px 10px" }}><span style={{ color: CLS_COLORS[r.cls], background: `${CLS_COLORS[r.cls]}18`, padding: "2px 8px", borderRadius: 4, fontSize: 11 }}>{r.cls}</span></td>
                <td style={{ padding: "6px 10px", color: sensitive ? "#ff4466" : "#ccd", fontWeight: sensitive ? 700 : 400, fontFamily: "monospace" }}>{r.pka.toFixed(2)}</td>
                {hasPropka && <td style={{ padding: "6px 10px", color: "#667788", fontFamily: "monospace" }}>{r.stdPka != null ? r.stdPka.toFixed(2) : "—"}</td>}
                {hasPropka && <td style={{ padding: "6px 10px", color: r.pkaShift != null && Math.abs(r.pkaShift) > 1 ? "#ff4466" : r.pkaShift != null && Math.abs(r.pkaShift) > 0.5 ? "#ffd43b" : "#667788", fontWeight: r.pkaShift != null && Math.abs(r.pkaShift) > 1 ? 700 : 400, fontFamily: "monospace" }}>{r.pkaShift != null ? `${r.pkaShift >= 0 ? "+" : ""}${r.pkaShift.toFixed(2)}` : "—"}</td>}
                <td style={{ padding: "6px 10px", color: sensitive ? "#ff4466" : "#8899aa", fontFamily: "monospace" }}>{delta >= 0 ? "+" : ""}{delta.toFixed(2)}</td>
                <td style={{ padding: "6px 10px" }}><span style={{ color: protonated ? "#4dabf7" : "#ff6b6b", fontSize: 11 }}>{state}</span></td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════
// Full Sequence Table
// ═══════════════════════════════════════════════════════════════
function FullSequenceTable({ sequence, titratable, ph }) {
  const titMap = {};
  for (const r of titratable) {
    titMap[r.pos] = r;
  }

  const SS_COLORS = {
    Acidic: "#ff6b6b", Basic: "#4dabf7", Thiol: "#ffd43b", Phenolic: "#da77f2"
  };

  return (
    <div style={{ maxHeight: 500, overflowY: "auto", borderRadius: 8, border: "1px solid #1e2a3a" }}>
      <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 12 }}>
        <thead>
          <tr style={{ position: "sticky", top: 0, background: "#0d1520", zIndex: 1 }}>
            {["Pos", "AA", "Amino Acid", "Titratable", "Class", "pKa", "State at pH"].map(h => (
              <th key={h} style={{ padding: "8px 10px", textAlign: "left", color: "#8899aa", borderBottom: "1px solid #1e2a3a", whiteSpace: "nowrap" }}>{h}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {sequence.split("").map((one, i) => {
            const pos = i + 1;
            const three = ONE_TO_THREE[one] || "???";
            const full = THREE_TO_FULL[three] || three;
            const tit = titMap[pos];
            const isTit = !!tit;
            const sensitive = tit && Math.abs(tit.pka - ph) < 1;
            const protonated = tit ? (ph < tit.pka) : null;

            return (
              <tr key={pos} style={{
                background: sensitive ? "#ff446610" : (pos % 2 === 0 ? "#0a1018" : "#0d1520"),
                borderBottom: "1px solid #111a25"
              }}>
                <td style={{ padding: "5px 10px", color: "#8899aa", fontFamily: "'JetBrains Mono', monospace", fontSize: 11 }}>{pos}</td>
                <td style={{ padding: "5px 10px", color: isTit ? (SS_COLORS[tit.cls] || "#ccd") : "#ccd", fontWeight: isTit ? 700 : 400, fontFamily: "'JetBrains Mono', monospace", fontSize: 13 }}>{one}</td>
                <td style={{ padding: "5px 10px", color: "#ccd" }}>{full}</td>
                <td style={{ padding: "5px 10px" }}>
                  {isTit
                    ? <span style={{ color: "#00cc88", fontSize: 11 }}>✓</span>
                    : <span style={{ color: "#334", fontSize: 11 }}>—</span>
                  }
                </td>
                <td style={{ padding: "5px 10px" }}>
                  {tit
                    ? <span style={{ color: SS_COLORS[tit.cls], background: `${SS_COLORS[tit.cls]}18`, padding: "2px 8px", borderRadius: 4, fontSize: 11 }}>{tit.cls}</span>
                    : <span style={{ color: "#334" }}>—</span>
                  }
                </td>
                <td style={{ padding: "5px 10px", fontFamily: "'JetBrains Mono', monospace", color: sensitive ? "#ff4466" : (isTit ? "#ccd" : "#334"), fontWeight: sensitive ? 700 : 400 }}>
                  {tit ? tit.pka.toFixed(2) : "—"}
                </td>
                <td style={{ padding: "5px 10px" }}>
                  {tit
                    ? <span style={{ color: protonated ? "#4dabf7" : "#ff6b6b", fontSize: 11 }}>{protonated ? "Protonated" : "Deprotonated"}</span>
                    : <span style={{ color: "#334" }}>—</span>
                  }
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════
// Molstar Viewer Component
// ═══════════════════════════════════════════════════════════════
const MolstarViewer = forwardRef(function MolstarViewer({ pdbId, customUrl, customData, customFormat }, ref) {
  const containerRef = useRef();
  const viewerRef = useRef(null);
  const [ready, setReady] = useState(false);

  // Expose viewer instance to parent
  useImperativeHandle(ref, () => ({
    getViewer: () => viewerRef.current,
    isReady: () => ready,
  }));

  useEffect(() => {
    if (!pdbId && !customUrl && !customData) return;

    const loadPdbeMolstar = () => {
      return new Promise((resolve) => {
        if (window.PDBeMolstarPlugin) { resolve(); return; }

        const css = document.createElement("link");
        css.rel = "stylesheet";
        css.href = "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.3.0/build/pdbe-molstar.css";
        document.head.appendChild(css);

        const js = document.createElement("script");
        js.src = "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.3.0/build/pdbe-molstar-plugin.js";
        js.onload = () => resolve();
        document.head.appendChild(js);
      });
    };

    loadPdbeMolstar().then(() => {
      if (viewerRef.current) {
        try { viewerRef.current.clear(); } catch (e) {}
      }

      setReady(false);
      const viewer = new PDBeMolstarPlugin();
      viewerRef.current = viewer;

      const options = {
        hideControls: false,
        hideCanvasControls: ["animation", "expand"],
        pdbeUrl: "https://www.ebi.ac.uk/pdbe/",
        bgColor: { r: 14, g: 17, b: 23 },
        highlightColor: { r: 255, g: 170, b: 51 },
        selectColor: { r: 0, g: 255, b: 136 },
        lighting: "metallic",
        sequencePanel: true,
        landscape: false,
        reactive: false,
        expanded: false,
      };

      if (pdbId) {
        options.moleculeId = pdbId.toLowerCase();
      } else if (customUrl) {
        options.customData = {
          url: customUrl,
          format: customFormat || "cif",
        };
      } else if (customData) {
        // For uploaded files, create a blob URL
        const blob = new Blob([customData], { type: "text/plain" });
        const blobUrl = URL.createObjectURL(blob);
        options.customData = {
          url: blobUrl,
          format: customFormat || "cif",
        };
      }

      viewer.render(containerRef.current, options);
      viewer.events.loadComplete.subscribe(() => setReady(true));
    });

    return () => {
      if (viewerRef.current) {
        try { viewerRef.current.clear(); } catch(e) {}
      }
    };
  }, [pdbId, customUrl, customData, customFormat]);

  return (
    <div style={{ position: "relative", width: "100%", height: 560, borderRadius: 8, overflow: "hidden", border: "1px solid #1e2a3a" }}>
      <div ref={containerRef} style={{ width: "100%", height: "100%" }} />
      {!ready && (
        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", background: "#0a1018ee", zIndex: 10 }}>
          <div style={{ textAlign: "center" }}>
            <div style={{ width: 32, height: 32, border: "3px solid #334", borderTopColor: "#4488ff", borderRadius: "50%", animation: "spin 1s linear infinite", margin: "0 auto 12px" }} />
            <div style={{ color: "#8899aa" }}>Loading structure...</div>
          </div>
        </div>
      )}
    </div>
  );
});

// ═══════════════════════════════════════════════════════════════
// First-Launch Onboarding — BRENDA Account Setup
// ═══════════════════════════════════════════════════════════════
function BrendaOnboarding({ onComplete, onSkip }) {
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [showPassword, setShowPassword] = useState(false);
  const [saving, setSaving] = useState(false);
  const [message, setMessage] = useState(null);
  const [step, setStep] = useState("intro"); // "intro" | "register" | "credentials"

  const handleSave = async () => {
    if (!email.trim() || !password.trim()) {
      setMessage({ type: "error", text: "Both email and password are required." });
      return;
    }
    setSaving(true);
    setMessage(null);
    const result = await saveBrendaCredentials(email.trim(), password);
    setSaving(false);
    if (result?.configured) {
      setMessage({ type: "success", text: "Credentials saved successfully!" });
      setTimeout(() => onComplete(), 800);
    } else {
      setMessage({ type: "error", text: result?.message || "Failed to save. Make sure the ConformSeek backend server is running." });
    }
  };

  const inputSt = { width: "100%", padding: "10px 12px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, color: "#d0d8e8", fontSize: 14, fontFamily: "'JetBrains Mono', monospace", outline: "none" };

  return (
    <div style={{ minHeight: "100vh", background: "#080c14", color: "#d0d8e8", fontFamily: "'DM Sans', 'Segoe UI', sans-serif", display: "flex", alignItems: "center", justifyContent: "center" }}>
      <div style={{ width: 520, padding: "0 20px" }}>

        {/* Logo */}
        <div style={{ textAlign: "center", marginBottom: 28 }}>
          <div style={{ width: 56, height: 56, borderRadius: 14, background: "linear-gradient(135deg, #4488ff, #00cc88)", display: "inline-flex", alignItems: "center", justifyContent: "center", fontSize: 28, marginBottom: 14 }}>🔬</div>
          <h1 style={{ fontSize: 28, fontWeight: 700, letterSpacing: -0.8, marginBottom: 4 }}>Welcome to ConformSeek</h1>
          <p style={{ fontSize: 13, color: "#556677" }}>pH-Dependent Protein Conformational Analysis</p>
        </div>

        {/* ── Step 1: Intro — Why BRENDA matters ── */}
        {step === "intro" && (
          <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 12, padding: "24px 24px 20px", marginBottom: 16 }}>
            <div style={{ fontSize: 15, fontWeight: 600, color: "#d0d8e8", marginBottom: 12 }}>Set up BRENDA for accurate enzyme pH data</div>
            <div style={{ fontSize: 13, color: "#8899aa", lineHeight: 1.7, marginBottom: 16 }}>
              <p style={{ marginBottom: 10 }}>
                ConformSeek uses the <span style={{ color: "#4488ff", fontWeight: 600 }}>BRENDA Enzyme Database</span> to
                look up experimentally determined <span style={{ color: "#00cc88", fontWeight: 600 }}>pH optima</span> for
                enzymes. This data comes from thousands of published studies and tells you at what pH each enzyme functions best.
              </p>
              <p style={{ marginBottom: 10 }}>
                With BRENDA data, ConformSeek can show you the <span style={{ color: "#00cc88" }}>functional pH</span> line
                on the titration curve and calculate how far a crystal structure's pH is from the enzyme's natural
                operating conditions — critical for understanding pH-driven conformational changes.
              </p>
              <p>
                BRENDA requires a <span style={{ color: "#ffd43b" }}>free account</span> to access their API.
                Your credentials are stored only on your computer and are sent directly to BRENDA's servers — never
                to any third party.
              </p>
            </div>

            <div style={{ display: "flex", gap: 10 }}>
              <button onClick={() => setStep("register")} style={{
                flex: 1, padding: "12px 0", background: "linear-gradient(135deg, #4488ff, #2266dd)", border: "none",
                borderRadius: 8, color: "#fff", fontSize: 14, fontWeight: 600, cursor: "pointer",
              }}>
                Set Up BRENDA Account
              </button>
              <button onClick={onSkip} style={{
                padding: "12px 20px", background: "transparent", border: "1px solid #1e2a3a",
                borderRadius: 8, color: "#556677", fontSize: 13, cursor: "pointer",
              }}>
                Skip for now
              </button>
            </div>

            <div style={{ fontSize: 11, color: "#445566", marginTop: 12, textAlign: "center", lineHeight: 1.5 }}>
              Without BRENDA, ConformSeek will still work — it will fall back to subcellular localization pH
              estimates from UniProt, which are less precise for enzymes.
            </div>
          </div>
        )}

        {/* ── Step 2: Registration instructions ── */}
        {step === "register" && (
          <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 12, padding: "24px 24px 20px", marginBottom: 16 }}>
            <div style={{ fontSize: 15, fontWeight: 600, color: "#d0d8e8", marginBottom: 12 }}>
              <span style={{ color: "#556677", fontSize: 11, textTransform: "uppercase", letterSpacing: "0.06em", display: "block", marginBottom: 4, fontWeight: 600 }}>Step 1 of 2</span>
              Register for a free BRENDA account
            </div>

            <div style={{ background: "#080c14", border: "1px solid #1e2a3a", borderRadius: 8, padding: "16px", marginBottom: 16 }}>
              <ol style={{ fontSize: 13, color: "#8899aa", lineHeight: 1.8, paddingLeft: 20, margin: 0 }}>
                <li style={{ marginBottom: 6 }}>
                  Visit the BRENDA registration page:
                  <br />
                  <a href="https://www.brenda-enzymes.org/register.php" target="_blank" rel="noopener noreferrer"
                    style={{ color: "#4488ff", textDecoration: "none", fontFamily: "'JetBrains Mono', monospace", fontSize: 12 }}>
                    brenda-enzymes.org/register.php ↗
                  </a>
                </li>
                <li style={{ marginBottom: 6 }}>
                  Enter your <span style={{ color: "#d0d8e8" }}>email address</span> and choose a <span style={{ color: "#d0d8e8" }}>password</span>
                </li>
                <li style={{ marginBottom: 6 }}>
                  Check your email for a <span style={{ color: "#ffd43b" }}>confirmation link</span> and activate your account
                </li>
                <li>
                  Come back here and enter your credentials below
                </li>
              </ol>
            </div>

            <div style={{ display: "flex", gap: 10 }}>
              <a href="https://www.brenda-enzymes.org/register.php" target="_blank" rel="noopener noreferrer"
                style={{
                  flex: 1, padding: "12px 0", background: "linear-gradient(135deg, #4488ff, #2266dd)", border: "none",
                  borderRadius: 8, color: "#fff", fontSize: 14, fontWeight: 600, cursor: "pointer", textAlign: "center",
                  textDecoration: "none", display: "block",
                }}>
                Open BRENDA Registration ↗
              </a>
              <button onClick={() => setStep("credentials")} style={{
                flex: 1, padding: "12px 0", background: "#00cc8820", border: "1px solid #00cc8840",
                borderRadius: 8, color: "#00cc88", fontSize: 14, fontWeight: 600, cursor: "pointer",
              }}>
                I have an account →
              </button>
            </div>

            <button onClick={() => setStep("intro")} style={{
              marginTop: 10, width: "100%", padding: "8px 0", background: "transparent", border: "none",
              color: "#556677", fontSize: 12, cursor: "pointer",
            }}>
              ← Back
            </button>
          </div>
        )}

        {/* ── Step 3: Enter credentials ── */}
        {step === "credentials" && (
          <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 12, padding: "24px 24px 20px", marginBottom: 16 }}>
            <div style={{ fontSize: 15, fontWeight: 600, color: "#d0d8e8", marginBottom: 12 }}>
              <span style={{ color: "#556677", fontSize: 11, textTransform: "uppercase", letterSpacing: "0.06em", display: "block", marginBottom: 4, fontWeight: 600 }}>Step 2 of 2</span>
              Enter your BRENDA credentials
            </div>

            <div style={{ display: "flex", flexDirection: "column", gap: 12, marginBottom: 16 }}>
              <div>
                <label style={{ fontSize: 12, color: "#8899aa", display: "block", marginBottom: 4 }}>Email address</label>
                <input type="email" value={email} onChange={e => setEmail(e.target.value)}
                  placeholder="you@example.edu" autoFocus style={inputSt} />
              </div>
              <div>
                <label style={{ fontSize: 12, color: "#8899aa", display: "block", marginBottom: 4 }}>Password</label>
                <div style={{ position: "relative" }}>
                  <input type={showPassword ? "text" : "password"} value={password} onChange={e => setPassword(e.target.value)}
                    placeholder="Your BRENDA password"
                    onKeyDown={e => { if (e.key === "Enter") handleSave(); }}
                    style={{ ...inputSt, paddingRight: 40 }} />
                  <button type="button" onClick={() => setShowPassword(!showPassword)}
                    style={{ position: "absolute", right: 8, top: "50%", transform: "translateY(-50%)", background: "none", border: "none", color: "#556677", fontSize: 16, cursor: "pointer" }}>
                    {showPassword ? "◉" : "◎"}
                  </button>
                </div>
              </div>
              <div style={{ fontSize: 11, color: "#445566", lineHeight: 1.4, padding: "8px 10px", background: "#080c14", borderRadius: 6, border: "1px solid #1a2030" }}>
                🔒 Your password is hashed (SHA-256) before storage — the plain text is never saved. Credentials
                are stored in <span style={{ fontFamily: "'JetBrains Mono', monospace", fontSize: 10 }}>~/.conformseek/</span> on your machine only.
              </div>
            </div>

            {message && (
              <div style={{
                marginBottom: 12, padding: "8px 14px", borderRadius: 6, fontSize: 12,
                background: message.type === "success" ? "#00cc8815" : "#ff446615",
                border: `1px solid ${message.type === "success" ? "#00cc8840" : "#ff446640"}`,
                color: message.type === "success" ? "#00cc88" : "#ff4466",
              }}>
                {message.text}
              </div>
            )}

            <div style={{ display: "flex", gap: 10 }}>
              <button onClick={handleSave} disabled={saving} style={{
                flex: 1, padding: "12px 0", background: "linear-gradient(135deg, #00cc88, #00aa66)", border: "none",
                borderRadius: 8, color: "#fff", fontSize: 14, fontWeight: 600,
                cursor: saving ? "wait" : "pointer", opacity: saving ? 0.6 : 1,
              }}>
                {saving ? "Saving..." : "Save & Continue"}
              </button>
              <button onClick={onSkip} style={{
                padding: "12px 20px", background: "transparent", border: "1px solid #1e2a3a",
                borderRadius: 8, color: "#556677", fontSize: 13, cursor: "pointer",
              }}>
                Skip
              </button>
            </div>

            <button onClick={() => setStep("register")} style={{
              marginTop: 10, width: "100%", padding: "8px 0", background: "transparent", border: "none",
              color: "#556677", fontSize: 12, cursor: "pointer",
            }}>
              ← I need to register first
            </button>
          </div>
        )}

        {/* Skip explanation footer (always visible) */}
        {step !== "intro" && (
          <div style={{ textAlign: "center", fontSize: 11, color: "#334455", lineHeight: 1.5 }}>
            You can always configure BRENDA credentials later via the ⚙ Settings menu.
          </div>
        )}
      </div>
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════
// Settings Panel (BRENDA credentials + cache management)
// ═══════════════════════════════════════════════════════════════
function SettingsPanel({ onClose, onCredentialsSaved }) {
  const [activeSection, setActiveSection] = useState("brenda");
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [saving, setSaving] = useState(false);
  const [credStatus, setCredStatus] = useState(null); // null = loading, object = loaded
  const [message, setMessage] = useState(null); // { type: "success"|"error", text }
  const [showPassword, setShowPassword] = useState(false);
  const [backendHealth, setBackendHealth] = useState(null);

  // Load credential status + backend health on mount
  useEffect(() => {
    fetchBrendaCredentialStatus().then(s => setCredStatus(s || { configured: false }));
    checkBackendHealth().then(h => setBackendHealth(h));
  }, []);

  const handleSave = async () => {
    if (!email.trim() || !password.trim()) {
      setMessage({ type: "error", text: "Both email and password are required." });
      return;
    }
    setSaving(true);
    setMessage(null);
    const result = await saveBrendaCredentials(email.trim(), password);
    setSaving(false);
    if (result?.configured) {
      setMessage({ type: "success", text: "Credentials saved! BRENDA pH lookups are now available." });
      setCredStatus(result);
      setEmail("");
      setPassword("");
      if (onCredentialsSaved) onCredentialsSaved();
    } else {
      setMessage({ type: "error", text: result?.message || "Failed to save credentials. Is the backend running?" });
    }
  };

  const handleDelete = async () => {
    setSaving(true);
    await deleteBrendaCredentials();
    setCredStatus({ configured: false });
    setMessage({ type: "success", text: "Credentials removed." });
    setSaving(false);
  };

  const handleClearCache = async () => {
    const result = await clearBrendaCache();
    if (result) {
      setMessage({ type: "success", text: "BRENDA query cache cleared." });
    }
  };

  const sectionBtn = (key, label) => ({
    padding: "8px 14px", fontSize: 12, fontWeight: activeSection === key ? 600 : 400,
    background: activeSection === key ? "#1e2a3a" : "transparent",
    border: "1px solid " + (activeSection === key ? "#2a3a4a" : "transparent"),
    borderRadius: 6, color: activeSection === key ? "#d0d8e8" : "#556677",
    cursor: "pointer", textAlign: "left", transition: "all 0.15s",
  });

  return (
    <div style={{ position: "fixed", inset: 0, zIndex: 9999, display: "flex", alignItems: "center", justifyContent: "center" }}>
      {/* Backdrop */}
      <div onClick={onClose} style={{ position: "absolute", inset: 0, background: "#000000aa", backdropFilter: "blur(4px)" }} />

      {/* Modal */}
      <div style={{ position: "relative", width: 580, maxHeight: "80vh", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 12, overflow: "hidden", display: "flex", flexDirection: "column" }}>

        {/* Header */}
        <div style={{ padding: "16px 20px", borderBottom: "1px solid #1e2a3a", display: "flex", alignItems: "center", justifyContent: "space-between" }}>
          <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
            <span style={{ fontSize: 18 }}>⚙</span>
            <span style={{ fontSize: 16, fontWeight: 600, color: "#d0d8e8" }}>Settings</span>
          </div>
          <button onClick={onClose} style={{ width: 28, height: 28, background: "#1e2a3a", border: "none", borderRadius: 6, color: "#8899aa", fontSize: 16, cursor: "pointer", display: "flex", alignItems: "center", justifyContent: "center" }}>✕</button>
        </div>

        <div style={{ display: "flex", flex: 1, overflow: "hidden" }}>
          {/* Sidebar */}
          <div style={{ width: 160, borderRight: "1px solid #1e2a3a", padding: "12px 8px", display: "flex", flexDirection: "column", gap: 4 }}>
            <button onClick={() => setActiveSection("brenda")} style={sectionBtn("brenda", "BRENDA API")}>BRENDA API</button>
            <button onClick={() => setActiveSection("about")} style={sectionBtn("about", "About")}>About</button>
          </div>

          {/* Content */}
          <div style={{ flex: 1, padding: "16px 20px", overflowY: "auto" }}>

            {activeSection === "brenda" && (
              <div>
                <div style={{ fontSize: 14, fontWeight: 600, color: "#d0d8e8", marginBottom: 4 }}>BRENDA Enzyme Database</div>
                <div style={{ fontSize: 12, color: "#556677", marginBottom: 16, lineHeight: 1.5 }}>
                  ConformSeek uses BRENDA's SOAP API to look up experimental pH optima for enzymes.
                  A free BRENDA account is required — your credentials are stored locally and never sent to anyone except BRENDA.
                </div>

                {/* Current status */}
                <div style={{ background: "#080c14", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 14px", marginBottom: 16 }}>
                  <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 8, fontWeight: 600 }}>Status</div>
                  <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
                    <div style={{ fontSize: 12, display: "flex", alignItems: "center", gap: 8 }}>
                      <span style={{ width: 8, height: 8, borderRadius: "50%", background: backendHealth ? "#00cc88" : "#ff4466", flexShrink: 0 }} />
                      <span style={{ color: "#8899aa" }}>Backend server: </span>
                      <span style={{ color: backendHealth ? "#00cc88" : "#ff4466" }}>{backendHealth ? "Connected" : "Not running"}</span>
                    </div>
                    <div style={{ fontSize: 12, display: "flex", alignItems: "center", gap: 8 }}>
                      <span style={{ width: 8, height: 8, borderRadius: "50%", background: backendHealth?.brenda_soap_available ? "#00cc88" : "#ffd43b", flexShrink: 0 }} />
                      <span style={{ color: "#8899aa" }}>Zeep SOAP library: </span>
                      <span style={{ color: backendHealth?.brenda_soap_available ? "#00cc88" : "#ffd43b" }}>
                        {backendHealth?.brenda_soap_available ? "Installed" : "Not installed"}
                      </span>
                    </div>
                    <div style={{ fontSize: 12, display: "flex", alignItems: "center", gap: 8 }}>
                      <span style={{ width: 8, height: 8, borderRadius: "50%", background: credStatus?.configured ? "#00cc88" : "#ff4466", flexShrink: 0 }} />
                      <span style={{ color: "#8899aa" }}>Credentials: </span>
                      <span style={{ color: credStatus?.configured ? "#00cc88" : "#ff4466" }}>
                        {credStatus === null ? "Checking..." : credStatus.configured ? credStatus.email : "Not configured"}
                      </span>
                    </div>
                    {backendHealth?.brenda_cached_ec_count > 0 && (
                      <div style={{ fontSize: 12, display: "flex", alignItems: "center", gap: 8 }}>
                        <span style={{ width: 8, height: 8, borderRadius: "50%", background: "#4488ff", flexShrink: 0 }} />
                        <span style={{ color: "#8899aa" }}>Cached queries: </span>
                        <span style={{ color: "#4488ff" }}>{backendHealth.brenda_cached_ec_count} EC numbers</span>
                      </div>
                    )}
                  </div>
                </div>

                {/* Credential form */}
                <div style={{ background: "#080c14", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 14px", marginBottom: 16 }}>
                  <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 8, fontWeight: 600 }}>
                    {credStatus?.configured ? "Update Credentials" : "Set Credentials"}
                  </div>
                  <div style={{ display: "flex", flexDirection: "column", gap: 8 }}>
                    <div>
                      <label style={{ fontSize: 11, color: "#8899aa", display: "block", marginBottom: 4 }}>Email</label>
                      <input
                        type="email" value={email} onChange={e => setEmail(e.target.value)}
                        placeholder={credStatus?.configured ? credStatus.email : "you@example.edu"}
                        style={{ width: "100%", padding: "8px 10px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 6, color: "#d0d8e8", fontSize: 13, fontFamily: "'JetBrains Mono', monospace", outline: "none" }}
                      />
                    </div>
                    <div>
                      <label style={{ fontSize: 11, color: "#8899aa", display: "block", marginBottom: 4 }}>Password</label>
                      <div style={{ position: "relative" }}>
                        <input
                          type={showPassword ? "text" : "password"} value={password} onChange={e => setPassword(e.target.value)}
                          placeholder="Your BRENDA password"
                          onKeyDown={e => { if (e.key === "Enter") handleSave(); }}
                          style={{ width: "100%", padding: "8px 36px 8px 10px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 6, color: "#d0d8e8", fontSize: 13, fontFamily: "'JetBrains Mono', monospace", outline: "none" }}
                        />
                        <button
                          type="button" onClick={() => setShowPassword(!showPassword)}
                          style={{ position: "absolute", right: 6, top: "50%", transform: "translateY(-50%)", background: "none", border: "none", color: "#556677", fontSize: 14, cursor: "pointer", padding: "2px 4px" }}
                        >{showPassword ? "◉" : "◎"}</button>
                      </div>
                    </div>
                    <div style={{ display: "flex", gap: 8, marginTop: 4 }}>
                      <button onClick={handleSave} disabled={saving}
                        style={{ flex: 1, padding: "8px 0", background: "linear-gradient(135deg, #4488ff, #2266dd)", border: "none", borderRadius: 6, color: "#fff", fontSize: 12, fontWeight: 600, cursor: saving ? "wait" : "pointer", opacity: saving ? 0.6 : 1 }}>
                        {saving ? "Saving..." : credStatus?.configured ? "Update Credentials" : "Save Credentials"}
                      </button>
                      {credStatus?.configured && (
                        <button onClick={handleDelete} disabled={saving}
                          style={{ padding: "8px 14px", background: "#ff446620", border: "1px solid #ff446640", borderRadius: 6, color: "#ff4466", fontSize: 12, fontWeight: 600, cursor: "pointer" }}>
                          Remove
                        </button>
                      )}
                    </div>
                  </div>
                  <div style={{ fontSize: 11, color: "#445566", marginTop: 8, lineHeight: 1.4 }}>
                    Don't have an account?{" "}
                    <a href="https://www.brenda-enzymes.org/register.php" target="_blank" rel="noopener noreferrer"
                      style={{ color: "#4488ff", textDecoration: "none" }}>
                      Register free at brenda-enzymes.org →
                    </a>
                  </div>
                </div>

                {/* Cache management */}
                {backendHealth?.brenda_cached_ec_count > 0 && (
                  <div style={{ background: "#080c14", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 14px" }}>
                    <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 8, fontWeight: 600 }}>Cache</div>
                    <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
                      <span style={{ fontSize: 12, color: "#8899aa" }}>{backendHealth.brenda_cached_ec_count} EC numbers cached locally</span>
                      <button onClick={handleClearCache} style={{ marginLeft: "auto", padding: "5px 12px", background: "#1e2a3a", border: "1px solid #2a3a4a", borderRadius: 4, color: "#8899aa", fontSize: 11, cursor: "pointer" }}>
                        Clear Cache
                      </button>
                    </div>
                  </div>
                )}

                {/* Status message */}
                {message && (
                  <div style={{
                    marginTop: 12, padding: "8px 14px", borderRadius: 6, fontSize: 12,
                    background: message.type === "success" ? "#00cc8815" : "#ff446615",
                    border: `1px solid ${message.type === "success" ? "#00cc8840" : "#ff446640"}`,
                    color: message.type === "success" ? "#00cc88" : "#ff4466",
                  }}>
                    {message.text}
                  </div>
                )}
              </div>
            )}

            {activeSection === "about" && (
              <div>
                <div style={{ fontSize: 14, fontWeight: 600, color: "#d0d8e8", marginBottom: 4 }}>ConformSeek</div>
                <div style={{ fontSize: 12, color: "#556677", lineHeight: 1.6 }}>
                  <p style={{ marginBottom: 12 }}>pH-Dependent Protein Conformational Analysis</p>
                  <p style={{ marginBottom: 8 }}>ConformSeek combines PROPKA pKa predictions with experimental enzyme pH data from BRENDA and subcellular localization data from UniProt to analyze how protein charge and conformation change with pH.</p>
                  <p style={{ marginBottom: 12 }}>Version 2.0 — BRENDA SOAP API (Bring Your Own Key)</p>
                  <div style={{ fontSize: 11, color: "#445566", borderTop: "1px solid #1e2a3a", paddingTop: 10, marginTop: 10 }}>
                    <div style={{ marginBottom: 4 }}>Author: Calvin (ConformSeek project)</div>
                    <div style={{ marginBottom: 4 }}>Backend: FastAPI + PROPKA3 + Zeep SOAP</div>
                    <div style={{ marginBottom: 4 }}>Frontend: React + D3.js + PDBe Molstar</div>
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════
// BRENDA Credentials Prompt (inline banner)
// ═══════════════════════════════════════════════════════════════
function BrendaCredentialsBanner({ onOpenSettings }) {
  return (
    <div style={{
      background: "#ffd43b10", border: "1px solid #ffd43b30", borderRadius: 8,
      padding: "10px 16px", marginBottom: 14, display: "flex", alignItems: "center", gap: 12,
    }}>
      <span style={{ fontSize: 18, flexShrink: 0 }}>🔑</span>
      <div style={{ flex: 1 }}>
        <div style={{ fontSize: 12, color: "#ffd43b", fontWeight: 600, marginBottom: 2 }}>BRENDA credentials required</div>
        <div style={{ fontSize: 11, color: "#8899aa", lineHeight: 1.4 }}>
          This enzyme has an EC number but pH optimum data couldn't be loaded.
          Set your free BRENDA account credentials in Settings to enable pH lookups.
        </div>
      </div>
      <button onClick={onOpenSettings} style={{
        padding: "6px 14px", background: "#ffd43b20", border: "1px solid #ffd43b40",
        borderRadius: 6, color: "#ffd43b", fontSize: 11, fontWeight: 600, cursor: "pointer", whiteSpace: "nowrap",
      }}>
        Open Settings
      </button>
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════
// Main App
// ═══════════════════════════════════════════════════════════════
export default function ConformSeekApp() {
  const [pdbId, setPdbId] = useState(null);
  const [inputId, setInputId] = useState("");
  const [inputMode, setInputMode] = useState("pdb");
  const [ph, setPh] = useState(7.4);
  const [info, setInfo] = useState(null);
  const [titratable, setTitratable] = useState([]);
  const [curve, setCurve] = useState([]);
  const [pI, setPI] = useState(null);
  const [loading, setLoading] = useState(false);
  const [activeTab, setActiveTab] = useState("titration");
  const [error, setError] = useState(null);
  const [uniprotResults, setUniprotResults] = useState(null);
  const [customUrl, setCustomUrl] = useState(null);
  const [customData, setCustomData] = useState(null);
  const [customFormat, setCustomFormat] = useState(null);
  const [uploadedFileName, setUploadedFileName] = useState(null);
  const [filePdbId, setFilePdbId] = useState(null);
  const [pkaSource, setPkaSource] = useState(null); // "propka" | "standard_fallback" | "approximation"
  const [brendaData, setBrendaData] = useState(null); // BRENDA pH optimum result
  const [subcellularPH, setSubcellularPH] = useState(null); // subcellular localization pH
  const [showSettings, setShowSettings] = useState(false); // Settings panel visibility
  const [brendaNeedsCredentials, setBrendaNeedsCredentials] = useState(false); // Show credential banner
  const [showOnboarding, setShowOnboarding] = useState(false); // First-launch onboarding
  const [onboardingChecked, setOnboardingChecked] = useState(false); // Whether we've checked yet
  const fileInputRef = useRef(null);
  const molstarRef = useRef(null);

  // ── First-launch onboarding check ──
  useEffect(() => {
    const checkFirstLaunch = async () => {
      const dismissed = localStorage.getItem("conformseek_onboarding_done");
      if (dismissed) {
        setOnboardingChecked(true);
        return;
      }
      // Check if credentials are already configured (maybe set via CLI/env)
      const status = await fetchBrendaCredentialStatus();
      if (status?.configured) {
        localStorage.setItem("conformseek_onboarding_done", "true");
        setOnboardingChecked(true);
        return;
      }
      // No credentials and never dismissed → show onboarding
      setShowOnboarding(true);
      setOnboardingChecked(true);
    };
    checkFirstLaunch();
  }, []);

  // ── Visualization layers ──
  const [layers, setLayers] = useState({
    ph_sensitive: true,
    titratable: false,
    protonation: false,
    pka_heatmap: false,
  });
  const toggleLayer = (key) => setLayers(prev => ({ ...prev, [key]: !prev[key] }));

  // ── Sequence range selection ──
  const [rangeStart, setRangeStart] = useState("");
  const [rangeEnd, setRangeEnd] = useState("");

  const loadProtein = useCallback(async (id) => {
    if (!id) return;
    setLoading(true);
    setError(null);
    setPkaSource(null);
    setBrendaData(null);
    setSubcellularPH(null);
    const data = await fetchPDBInfo(id);
    setInfo(data);

    // Try PROPKA backend first, fall back to client-side approximation
    let tit = null;
    const propkaData = await fetchPropkaFromBackend(id);
    console.log("PROPKA result:", propkaData ? `${propkaData.titratable.length} residues (${propkaData.method})` : "unavailable");
    if (propkaData && propkaData.titratable.length > 0) {
      tit = propkaToTitratable(propkaData);
      setPkaSource(propkaData.method); // "propka" or "standard_fallback"
      setPI(propkaData.pI);
    } else if (data.sequence) {
      tit = buildTitratableFromSequence(data.sequence);
      setPkaSource("approximation");
      setPI(computePI(tit));
    }

    if (tit) {
      setTitratable(tit);
      setCurve(computeTitrationCurve(tit));
    }

    // Fetch BRENDA pH data if we have an EC number
    // Then fall back to subcellular localization pH for non-enzymes
    console.log("Functional pH lookup:", { ecNumber: data.ecNumber, uniprotAcc: data.uniprotAcc });
    setBrendaNeedsCredentials(false);
    if (data.ecNumber) {
      fetchBrendaPH(data.ecNumber).then(bd => {
        console.log("BRENDA result:", bd ? (bd._needsCredentials ? "needs credentials" : `EC ${bd.ec_number}, pH ${bd.ph_median}`) : "no data");
        if (bd && bd._needsCredentials) {
          setBrendaNeedsCredentials(true);
          // Still try subcellular as fallback
          if (data.uniprotAcc) {
            fetchSubcellularPH(data.uniprotAcc).then(sc => {
              console.log("Subcellular pH result:", sc || "no data");
              if (sc) setSubcellularPH(sc);
            });
          }
        } else if (bd) {
          setBrendaData(bd);
        } else if (data.uniprotAcc) {
          fetchSubcellularPH(data.uniprotAcc).then(sc => {
            console.log("Subcellular pH result:", sc || "no data");
            if (sc) setSubcellularPH(sc);
          });
        }
      });
    } else if (data.uniprotAcc) {
      // No EC number — use subcellular localization (tier 2)
      fetchSubcellularPH(data.uniprotAcc).then(sc => {
        console.log("Subcellular pH result:", sc || "no data");
        if (sc) setSubcellularPH(sc);
      });
    }

    setLoading(false);
  }, []);

  useEffect(() => { if (pdbId) loadProtein(pdbId); }, [pdbId, loadProtein]);

  // ── Apply visual overlays to Mol* when layers/pH change ──
  useEffect(() => {
    const viewer = molstarRef.current?.getViewer();
    if (!viewer || !molstarRef.current?.isReady()) return;

    const applyVisuals = async () => {
      try {
        await viewer.visual.clearSelection(undefined, { keepColors: false, keepRepresentations: false });
      } catch (e) {}

      const data = [];

      if (layers.ph_sensitive) {
        for (const r of titratable) {
          if (Math.abs(r.pka - ph) < 1.0) {
            data.push({
              struct_asym_id: "A",
              start_residue_number: r.pos,
              end_residue_number: r.pos,
              color: "#ff4466",
              sideChain: true,
            });
          }
        }
      }

      if (layers.titratable) {
        for (const r of titratable) {
          if (layers.ph_sensitive && Math.abs(r.pka - ph) < 1.0) continue;
          const clsColor = { Acidic: "#ff6b6b", Basic: "#4dabf7", Thiol: "#ffd43b", Phenolic: "#da77f2" };
          data.push({
            struct_asym_id: "A",
            start_residue_number: r.pos,
            end_residue_number: r.pos,
            color: clsColor[r.cls] || "#88aacc",
            sideChain: true,
          });
        }
      }

      if (layers.protonation) {
        for (const r of titratable) {
          if (data.some(d => d.start_residue_number === r.pos)) continue;
          const protonated = ph < r.pka;
          data.push({
            struct_asym_id: "A",
            start_residue_number: r.pos,
            end_residue_number: r.pos,
            color: protonated ? "#4dabf7" : "#ff6b6b",
            sideChain: true,
          });
        }
      }

      if (layers.pka_heatmap) {
        for (const r of titratable) {
          if (data.some(d => d.start_residue_number === r.pos)) continue;
          const dist = Math.abs(r.pka - ph);
          const t = Math.min(dist / 7, 1);
          const red = Math.round(255 * (1 - t));
          const blue = Math.round(255 * t);
          const hex = `#${red.toString(16).padStart(2, "0")}44${blue.toString(16).padStart(2, "0")}`;
          data.push({
            struct_asym_id: "A",
            start_residue_number: r.pos,
            end_residue_number: r.pos,
            color: hex,
            sideChain: true,
          });
        }
      }

      if (data.length > 0) {
        try {
          // PDBe Molstar can choke on large data arrays — batch if needed
          const BATCH = 150;
          if (data.length <= BATCH) {
            await viewer.visual.select({ data, nonSelectedColor: undefined });
          } else {
            for (let i = 0; i < data.length; i += BATCH) {
              const chunk = data.slice(i, i + BATCH);
              await viewer.visual.select({ data: chunk, nonSelectedColor: undefined });
            }
          }
        } catch (e) { console.warn("Mol* select failed:", e); }
      }
    };

    // Small delay to ensure Mol* is ready
    const timer = setTimeout(applyVisuals, 300);
    return () => clearTimeout(timer);
  }, [layers, ph, titratable]);

  // ── Apply sequence range focus ──
  const applyRangeFocus = () => {
    const viewer = molstarRef.current?.getViewer();
    if (!viewer || !molstarRef.current?.isReady()) return;
    const start = parseInt(rangeStart);
    const end = parseInt(rangeEnd);
    if (isNaN(start) || isNaN(end) || start > end) return;

    viewer.visual.focus([{
      struct_asym_id: "A",
      start_residue_number: start,
      end_residue_number: end,
    }]);

    viewer.visual.select({
      data: [{
        struct_asym_id: "A",
        start_residue_number: start,
        end_residue_number: end,
        color: "#00ff88",
        sideChain: true,
      }],
    });
  };

  const clearRange = () => {
    const viewer = molstarRef.current?.getViewer();
    if (viewer) {
      try { viewer.visual.clearSelection(); } catch (e) {}
      try { viewer.visual.reset({ camera: true }); } catch (e) {}
    }
    setRangeStart("");
    setRangeEnd("");
  };

  // ── Quick select functions ──
  const quickSelectTitratable = () => {
    const viewer = molstarRef.current?.getViewer();
    if (!viewer || !molstarRef.current?.isReady()) return;
    const clsColor = { Acidic: "#ff6b6b", Basic: "#4dabf7", Thiol: "#ffd43b", Phenolic: "#da77f2" };
    const data = titratable.map(r => ({
      struct_asym_id: "A",
      start_residue_number: r.pos,
      end_residue_number: r.pos,
      color: clsColor[r.cls] || "#88aacc",
      sideChain: true,
    }));
    viewer.visual.select({ data });
  };

  const quickSelectPhSensitive = () => {
    const viewer = molstarRef.current?.getViewer();
    if (!viewer || !molstarRef.current?.isReady()) return;
    const data = titratable.filter(r => Math.abs(r.pka - ph) < 1.0).map(r => ({
      struct_asym_id: "A",
      start_residue_number: r.pos,
      end_residue_number: r.pos,
      color: "#ff4466",
      sideChain: true,
      focus: false,
    }));
    if (data.length > 0) viewer.visual.select({ data });
  };

  const quickClearSelection = () => {
    const viewer = molstarRef.current?.getViewer();
    if (viewer) {
      try { viewer.visual.clearSelection(); } catch (e) {}
    }
  };

  const resetState = () => {
    setError(null); setUniprotResults(null);
    setCustomUrl(null); setCustomData(null); setCustomFormat(null);
    setUploadedFileName(null); setFilePdbId(null);
    setBrendaNeedsCredentials(false);
  };

  // When credentials are saved in settings, retry the BRENDA lookup if needed
  const handleCredentialsSaved = () => {
    setBrendaNeedsCredentials(false);
    if (info?.ecNumber && !brendaData) {
      fetchBrendaPH(info.ecNumber).then(bd => {
        if (bd && !bd._needsCredentials) setBrendaData(bd);
      });
    }
  };

  const handlePdbSubmit = (e) => {
    e.preventDefault();
    resetState();
    const id = inputId.trim().toUpperCase();
    if (id.length >= 4) setPdbId(id);
  };

  const handleUniprotSubmit = async (e) => {
    e.preventDefault();
    resetState();
    setLoading(true);
    setPdbId(null);
    const id = inputId.trim();

    const result = await lookupUniProt(id);

    if (result.type === "pdb_list") {
      if (result.pdbIds.length === 1) {
        setPdbId(result.pdbIds[0]);
        setLoading(false);
      } else {
        setUniprotResults(result);
        setLoading(false);
      }
    } else if (result.type === "alphafold") {
      setCustomUrl(result.url);
      setCustomFormat("cif");
      setInfo({
        title: `AlphaFold Predicted Structure — ${result.entryId}`,
        method: "AlphaFold Prediction", resolution: null, crystPh: null,
        crystTemp: null, crystMethod: null, crystDetails: null,
        collectionTemp: null, rFree: null, rWork: null,
        deposited: null, released: null, keywords: "Predicted structure",
        sequence: "", organism: null, expressionSystem: null,
      });
      try {
        const seqRes = await fetch(`https://rest.uniprot.org/uniprotkb/${id}.fasta`);
        if (seqRes.ok) {
          const fasta = await seqRes.text();
          const seq = fasta.split("\n").filter(l => !l.startsWith(">")).join("");
          const tit = buildTitratableFromSequence(seq);
          setTitratable(tit);
          setCurve(computeTitrationCurve(tit));
          setPI(computePI(tit));
          setPkaSource("approximation");
          setInfo(prev => ({ ...prev, sequence: seq }));
        }
      } catch (err) { console.error("UniProt FASTA fetch failed:", err); }
      setLoading(false);
    } else {
      setError(`No structures found for UniProt ID: ${id}`);
      setLoading(false);
    }
  };

// ═══════════════════════════════════════════════════════════════
// Parse metadata from mmCIF/PDB file text
// ═══════════════════════════════════════════════════════════════
function parseCifValue(text, category, field) {
  // Handles: _category.field 'value', _category.field value, and multi-line ;value;
  const key = `${category}.${field}`;
  // Single-line quoted: _key 'value' or _key "value"
  const m1 = text.match(new RegExp(`${key.replace('.', '\\.')}\\s+['"](.*?)['"]`, 's'));
  if (m1) return m1[1].trim();
  // Single-line unquoted: _key value
  const m2 = text.match(new RegExp(`${key.replace('.', '\\.')}\\s+(\\S+)`));
  if (m2 && m2[1] !== '?' && m2[1] !== '.') return m2[1].trim();
  // Multi-line semicolon delimited: _key\n;value\n;
  const m3 = text.match(new RegExp(`${key.replace('.', '\\.')}\\s*\\n;([^;]+);`, 's'));
  if (m3) return m3[1].trim();
  return null;
}

function parseMetadataFromFile(text, format) {
  const meta = {};

  if (format === "cif") {
    meta.pdbId = parseCifValue(text, '_entry', 'id');
    meta.title = parseCifValue(text, '_struct', 'title');
    meta.method = parseCifValue(text, '_exptl', 'method');
    meta.resolution = parseCifValue(text, '_refine', 'ls_d_res_high');
    meta.crystPh = parseCifValue(text, '_exptl_crystal_grow', 'pH');
    meta.crystTemp = parseCifValue(text, '_exptl_crystal_grow', 'temp');
    meta.crystMethod = parseCifValue(text, '_exptl_crystal_grow', 'method');
    meta.crystDetails = parseCifValue(text, '_exptl_crystal_grow', 'pdbx_details');
    meta.collectionTemp = parseCifValue(text, '_diffrn', 'ambient_temp');
    meta.rFree = parseCifValue(text, '_refine', 'ls_R_factor_R_free');
    meta.rWork = parseCifValue(text, '_refine', 'ls_R_factor_R_work');
    meta.keywords = parseCifValue(text, '_struct_keywords', 'pdbx_keywords');
    // Organism from _entity_src_gen or _pdbx_entity_src_syn
    meta.organism = parseCifValue(text, '_entity_src_gen', 'pdbx_gene_src_scientific_name')
                 || parseCifValue(text, '_entity_src_nat', 'pdbx_organism_scientific');
    meta.expressionSystem = parseCifValue(text, '_entity_src_gen', 'pdbx_host_org_scientific_name');
  } else {
    // PDB format
    const lines = text.split("\n");
    const titleParts = [];
    for (const line of lines) {
      if (line.startsWith("HEADER")) {
        meta.pdbId = line.substring(62, 66).trim();
        meta.keywords = line.substring(10, 50).trim();
      }
      if (line.startsWith("TITLE")) titleParts.push(line.substring(10).trim());
      if (line.startsWith("EXPDTA")) meta.method = line.substring(10).trim();
      if (line.startsWith("REMARK   2 RESOLUTION.")) {
        const rm = line.match(/([\d.]+)\s*ANGSTROMS/);
        if (rm) meta.resolution = rm[1];
      }
      if (line.startsWith("SOURCE")) {
        const orgMatch = line.match(/ORGANISM_SCIENTIFIC:\s*(.+?)(?:;|$)/);
        if (orgMatch) meta.organism = orgMatch[1].trim();
        const expMatch = line.match(/EXPRESSION_SYSTEM:\s*(.+?)(?:;|$)/);
        if (expMatch) meta.expressionSystem = expMatch[1].trim();
      }
      if (line.startsWith("ATOM")) break;
    }
    if (titleParts.length) meta.title = titleParts.join(" ");
  }

  return meta;
}

  const handleFileUpload = async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    resetState();
    setPdbId(null);

    const text = await file.text();
    const name = file.name.toLowerCase();
    const format = (name.endsWith(".cif") || name.endsWith(".mmcif")) ? "cif" : "pdb";

    setCustomData(text);
    setCustomFormat(format);
    setUploadedFileName(file.name);

    // Parse sequence
    const seq = parseSequenceFromFileText(text, format);
    const tit = buildTitratableFromSequence(seq);
    setTitratable(tit);
    setCurve(computeTitrationCurve(tit));
    setPI(computePI(tit));
    setPkaSource("approximation");

    // Parse metadata from the file itself
    const fileMeta = parseMetadataFromFile(text, format);

    const toFloat = (v) => { const n = parseFloat(v); return isNaN(n) ? null : n; };

    let infoObj = {
      title: fileMeta.title || file.name.replace(/\.(pdb|cif|mmcif|ent)$/i, ""),
      method: fileMeta.method || "Uploaded file",
      resolution: toFloat(fileMeta.resolution),
      crystPh: toFloat(fileMeta.crystPh),
      crystTemp: toFloat(fileMeta.crystTemp),
      crystMethod: fileMeta.crystMethod || null,
      crystDetails: fileMeta.crystDetails || null,
      collectionTemp: toFloat(fileMeta.collectionTemp),
      rFree: toFloat(fileMeta.rFree),
      rWork: toFloat(fileMeta.rWork),
      deposited: null, released: null,
      keywords: fileMeta.keywords || null,
      sequence: seq,
      organism: fileMeta.organism || null,
      expressionSystem: fileMeta.expressionSystem || null,
    };

    // If we found a PDB ID in the file, try to fetch full metadata from RCSB API
    if (fileMeta.pdbId && fileMeta.pdbId.length >= 4) {
      setFilePdbId(fileMeta.pdbId.toUpperCase());
      try {
        const apiData = await fetchPDBInfo(fileMeta.pdbId);
        if (apiData.title) {
          // Merge: prefer API data but keep file-parsed sequence
          infoObj = { ...apiData, sequence: seq.length > 0 ? seq : apiData.sequence };
        }
      } catch (err) {
        console.log("RCSB API fetch for uploaded file PDB ID failed, using file metadata:", err);
      }
    }

    setInfo(infoObj);
  };

  const handleSelectPdb = (id) => {
    setUniprotResults(null);
    setPdbId(id);
  };

  const hasStructure = pdbId || customUrl || customData;

  const nSensitive = titratable.filter(r => Math.abs(r.pka - ph) < 1).length;
  const nAcidic = titratable.filter(r => r.cls === "Acidic").length;
  const nBasic = titratable.filter(r => r.cls === "Basic").length;
  const charge = curve.length > 0 ? curve.reduce((p, c) => Math.abs(c.ph - ph) < Math.abs(p.ph - ph) ? c : p).charge : 0;

  // ── Shared CSS ──
  const globalCSS = `
    @import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap');
    @keyframes spin { to { transform: rotate(360deg); } }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    ::-webkit-scrollbar { width: 6px; height: 6px; }
    ::-webkit-scrollbar-track { background: #0a1018; }
    ::-webkit-scrollbar-thumb { background: #1e2a3a; border-radius: 3px; }
    input[type="range"] { -webkit-appearance: none; appearance: none; background: transparent; width: 100%; height: 20px; cursor: pointer; }
    input[type="range"]::-webkit-slider-runnable-track { height: 6px; background: linear-gradient(90deg, #ff4466, #ffcc00, #4488ff); border-radius: 3px; border: 1px solid #334; }
    input[type="range"]::-webkit-slider-thumb { -webkit-appearance: none; appearance: none; width: 20px; height: 20px; background: #fff; border-radius: 50%; margin-top: -8px; box-shadow: 0 0 8px #ff446688; cursor: pointer; }
    input[type="range"]::-moz-range-track { height: 6px; background: linear-gradient(90deg, #ff4466, #ffcc00, #4488ff); border-radius: 3px; border: 1px solid #334; }
    input[type="range"]::-moz-range-thumb { width: 20px; height: 20px; background: #fff; border-radius: 50%; border: none; box-shadow: 0 0 8px #ff446688; cursor: pointer; }
    input[type="range"]::-moz-range-progress { background: transparent; }
    .titlebar-drag { -webkit-app-region: drag; }
    .titlebar-no-drag { -webkit-app-region: no-drag; }
  `;

  // Shared input styles
  const inputStyle = { padding: "10px 14px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, color: "#d0d8e8", fontSize: 14, fontFamily: "'JetBrains Mono', monospace", textAlign: "center", outline: "none" };
  const btnStyle = { padding: "10px 24px", background: "linear-gradient(135deg, #4488ff, #2266dd)", border: "none", borderRadius: 8, color: "#fff", fontSize: 14, fontWeight: 600, cursor: "pointer" };
  const tabBtnStyle = (active) => ({ padding: "8px 18px", background: active ? "#4488ff15" : "transparent", border: `1px solid ${active ? "#4488ff" : "#1e2a3a"}`, borderRadius: 8, color: active ? "#4488ff" : "#556677", fontSize: 13, fontWeight: active ? 600 : 400, cursor: "pointer", transition: "all 0.15s" });

  // ═══ First-launch onboarding ═══
  if (showOnboarding) {
    return (
      <>
        <style>{globalCSS}</style>
        <BrendaOnboarding
          onComplete={() => {
            localStorage.setItem("conformseek_onboarding_done", "true");
            setShowOnboarding(false);
          }}
          onSkip={() => {
            localStorage.setItem("conformseek_onboarding_done", "true");
            setShowOnboarding(false);
          }}
        />
      </>
    );
  }

  // ═══ Landing page ═══
  if (!hasStructure && !uniprotResults) {
    return (
      <div style={{ minHeight: "100vh", background: "#080c14", color: "#d0d8e8", fontFamily: "'DM Sans', 'Segoe UI', sans-serif" }}>
        <style>{globalCSS}</style>
        {/* Draggable title bar region */}
        <div className="titlebar-drag" style={{ height: 38, width: "100%", position: "fixed", top: 0, left: 0, zIndex: 100 }} />
        <div className="titlebar-no-drag" style={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center", minHeight: "100vh", padding: 40 }}>
          <div style={{ width: 64, height: 64, borderRadius: 16, background: "linear-gradient(135deg, #4488ff, #00cc88)", display: "flex", alignItems: "center", justifyContent: "center", fontSize: 32, marginBottom: 20 }}>🔬</div>
          <h1 style={{ fontSize: 36, fontWeight: 700, letterSpacing: -1, marginBottom: 8 }}>ConformSeek</h1>
          <p style={{ fontSize: 14, color: "#556677", marginBottom: 32 }}>pH-Dependent Protein Conformational Analysis</p>

          {/* Mode tabs */}
          <div style={{ display: "flex", gap: 8, marginBottom: 20 }}>
            {[["pdb", "PDB ID"], ["uniprot", "UniProt ID"], ["upload", "Upload File"]].map(([mode, label]) => (
              <button key={mode} onClick={() => { setInputMode(mode); setError(null); }} style={tabBtnStyle(inputMode === mode)}>{label}</button>
            ))}
          </div>

          {/* PDB ID input */}
          {inputMode === "pdb" && (
            <form onSubmit={handlePdbSubmit} style={{ display: "flex", gap: 10, marginBottom: 12 }}>
              <input value={inputId} onChange={e => setInputId(e.target.value.toUpperCase())} placeholder="e.g. 3CRI, 9NXE" autoFocus style={{ ...inputStyle, width: 240 }} />
              <button type="submit" style={btnStyle} disabled={loading}>{loading ? "Loading..." : "Load"}</button>
            </form>
          )}

          {/* UniProt ID input */}
          {inputMode === "uniprot" && (
            <form onSubmit={handleUniprotSubmit} style={{ display: "flex", gap: 10, marginBottom: 12 }}>
              <input value={inputId} onChange={e => setInputId(e.target.value)} placeholder="e.g. P00760, Q9Y5Y9" autoFocus style={{ ...inputStyle, width: 240 }} />
              <button type="submit" style={btnStyle} disabled={loading}>{loading ? "Searching..." : "Search"}</button>
            </form>
          )}

          {/* File upload */}
          {inputMode === "upload" && (
            <div style={{ textAlign: "center", marginBottom: 12 }}>
              <input ref={fileInputRef} type="file" accept=".pdb,.cif,.mmcif,.ent" onChange={handleFileUpload} style={{ display: "none" }} />
              <button onClick={() => fileInputRef.current?.click()} style={{ ...btnStyle, padding: "14px 32px", background: "#0d1520", border: "1px solid #1e2a3a", color: "#d0d8e8" }}>
                📁 Choose PDB or mmCIF file
              </button>
              <p style={{ fontSize: 11, color: "#334455", marginTop: 8 }}>Supports .pdb, .cif, .mmcif, .ent</p>
            </div>
          )}

          {/* Hint text */}
          <p style={{ fontSize: 12, color: "#334455", marginTop: 4 }}>
            {inputMode === "pdb" && "Enter a 4-character PDB ID to load an experimental structure"}
            {inputMode === "uniprot" && "Enter a UniProt accession — maps to PDB structures or AlphaFold prediction"}
            {inputMode === "upload" && "Upload a local structure file for analysis"}
          </p>

          {/* Settings link */}
          <button onClick={() => setShowSettings(true)} style={{ marginTop: 20, padding: "8px 18px", background: "transparent", border: "1px solid #1e2a3a", borderRadius: 8, color: "#556677", fontSize: 12, cursor: "pointer", display: "flex", alignItems: "center", gap: 6, transition: "all 0.15s" }}
            onMouseEnter={e => { e.target.style.borderColor = "#2a3a4a"; e.target.style.color = "#8899aa"; }}
            onMouseLeave={e => { e.target.style.borderColor = "#1e2a3a"; e.target.style.color = "#556677"; }}
          >
            <span>⚙</span> Settings — Configure BRENDA API credentials
          </button>

          {/* Settings Modal (landing page) */}
          {showSettings && <SettingsPanel onClose={() => setShowSettings(false)} onCredentialsSaved={handleCredentialsSaved} />}

          {/* Error */}
          {error && (
            <div style={{ marginTop: 16, padding: "10px 20px", background: "#ff446620", border: "1px solid #ff446650", borderRadius: 8, color: "#ff6b6b", fontSize: 13 }}>{error}</div>
          )}

          {/* Loading spinner */}
          {loading && (
            <div style={{ marginTop: 20 }}>
              <div style={{ width: 28, height: 28, border: "3px solid #334", borderTopColor: "#4488ff", borderRadius: "50%", animation: "spin 1s linear infinite", margin: "0 auto" }} />
            </div>
          )}
        </div>
      </div>
    );
  }

  // ═══ UniProt results picker (multiple PDB hits) ═══
  if (uniprotResults && !hasStructure) {
    return (
      <div style={{ minHeight: "100vh", background: "#080c14", color: "#d0d8e8", fontFamily: "'DM Sans', 'Segoe UI', sans-serif" }}>
        <style>{globalCSS}</style>
        <div className="titlebar-drag" style={{ height: 38, width: "100%", position: "fixed", top: 0, left: 0, zIndex: 100 }} />
        <div className="titlebar-no-drag" style={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center", minHeight: "100vh", padding: 40 }}>
          <div style={{ width: 48, height: 48, borderRadius: 12, background: "linear-gradient(135deg, #4488ff, #00cc88)", display: "flex", alignItems: "center", justifyContent: "center", fontSize: 24, marginBottom: 16 }}>🔬</div>
          <h2 style={{ fontSize: 22, fontWeight: 600, marginBottom: 6 }}>Multiple PDB structures found</h2>
          <p style={{ fontSize: 13, color: "#556677", marginBottom: 20 }}>
            UniProt <span style={{ color: "#4488ff", fontFamily: "'JetBrains Mono', monospace" }}>{uniprotResults.uniprotId}</span> maps to {uniprotResults.pdbIds.length} PDB entries — select one:
          </p>
          <div style={{ display: "flex", flexWrap: "wrap", gap: 8, justifyContent: "center", maxWidth: 600, marginBottom: 20 }}>
            {uniprotResults.pdbIds.map(id => (
              <button key={id} onClick={() => handleSelectPdb(id)} style={{ padding: "10px 20px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, color: "#4488ff", fontSize: 14, fontFamily: "'JetBrains Mono', monospace", fontWeight: 600, cursor: "pointer", transition: "all 0.15s" }}
                onMouseEnter={e => { e.target.style.borderColor = "#4488ff"; e.target.style.background = "#4488ff15"; }}
                onMouseLeave={e => { e.target.style.borderColor = "#1e2a3a"; e.target.style.background = "#0d1520"; }}
              >{id}</button>
            ))}
          </div>
          <button onClick={() => { setUniprotResults(null); setPdbId(null); }} style={{ padding: "8px 20px", background: "transparent", border: "1px solid #334", borderRadius: 6, color: "#556677", fontSize: 12, cursor: "pointer" }}>← Back</button>
        </div>
      </div>
    );
  }

  return (
    <div style={{ minHeight: "100vh", background: "#080c14", color: "#d0d8e8", fontFamily: "'DM Sans', 'Segoe UI', sans-serif" }}>
      <style>{globalCSS}</style>

      {/* Settings Modal */}
      {showSettings && <SettingsPanel onClose={() => setShowSettings(false)} onCredentialsSaved={handleCredentialsSaved} />}
      <header className="titlebar-drag" style={{ padding: "12px 28px 12px 84px", borderBottom: "1px solid #1e2a3a", display: "flex", alignItems: "center", gap: 16, background: "#0a0e18" }}>
        <div className="titlebar-no-drag" onClick={() => { setPdbId(null); setCustomUrl(null); setCustomData(null); setInfo(null); setTitratable([]); setCurve([]); setPI(null); setPkaSource(null); setBrendaData(null); setSubcellularPH(null); setUniprotResults(null); setUploadedFileName(null); setBrendaNeedsCredentials(false); }} style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
          <div style={{ width: 28, height: 28, borderRadius: 6, background: "linear-gradient(135deg, #4488ff, #00cc88)", display: "flex", alignItems: "center", justifyContent: "center", fontSize: 15 }}>🔬</div>
          <span style={{ fontSize: 16, fontWeight: 700, letterSpacing: -0.5 }}>ConformSeek</span>
        </div>
        <form className="titlebar-no-drag" onSubmit={handlePdbSubmit} style={{ display: "flex", gap: 6, marginLeft: "auto" }}>
          <input value={inputId} onChange={e => setInputId(e.target.value)} placeholder="PDB or UniProt ID" style={{ width: 170, padding: "6px 10px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 6, color: "#d0d8e8", fontSize: 12, fontFamily: "'JetBrains Mono', monospace" }} />
          <button type="submit" style={{ padding: "6px 14px", background: "linear-gradient(135deg, #4488ff, #2266dd)", border: "none", borderRadius: 6, color: "#fff", fontSize: 12, fontWeight: 600, cursor: "pointer" }}>PDB</button>
          <button type="button" onClick={handleUniprotSubmit} style={{ padding: "6px 14px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 6, color: "#00cc88", fontSize: 12, fontWeight: 600, cursor: "pointer" }}>UniProt</button>
          <input ref={fileInputRef} type="file" accept=".pdb,.cif,.mmcif,.ent" onChange={handleFileUpload} style={{ display: "none" }} />
          <button type="button" onClick={() => fileInputRef.current?.click()} style={{ padding: "6px 14px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 6, color: "#ffd43b", fontSize: 12, fontWeight: 600, cursor: "pointer" }}>📁 Upload</button>
          <button type="button" onClick={() => setShowSettings(true)} style={{ padding: "6px 10px", background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 6, color: "#8899aa", fontSize: 14, cursor: "pointer", marginLeft: 4, transition: "all 0.15s" }}
            onMouseEnter={e => { e.target.style.color = "#d0d8e8"; e.target.style.borderColor = "#2a3a4a"; }}
            onMouseLeave={e => { e.target.style.color = "#8899aa"; e.target.style.borderColor = "#1e2a3a"; }}
            title="Settings"
          >⚙</button>
        </form>
      </header>

      <div style={{ maxWidth: 1280, margin: "0 auto", padding: "20px 24px" }}>
        {/* ── Title ── */}
        {info?.title && (
          <div style={{ marginBottom: 12 }}>
            <h1 style={{ fontSize: 22, fontWeight: 600, color: "#e8ecf4", lineHeight: 1.3, letterSpacing: -0.3 }}>{info.title}</h1>
            <div style={{ fontSize: 12, color: "#556677", marginTop: 4 }}>
              {pdbId && <>PDB: <span style={{ color: "#4488ff", fontFamily: "'JetBrains Mono', monospace" }}>{pdbId}</span></>}
              {uploadedFileName && <>File: <span style={{ color: "#ffd43b", fontFamily: "'JetBrains Mono', monospace" }}>{uploadedFileName}</span></>}
              {!pdbId && customUrl && !uploadedFileName && <>Source: <span style={{ color: "#00cc88" }}>AlphaFold</span></>}
              {info.organism && <> · {info.organism}</>}
              {info.keywords && <> · {info.keywords}</>}
            </div>
          </div>
        )}

        {/* ── Experimental Conditions ── */}
        {info && (
          <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 16px", marginBottom: 14, display: "flex", flexWrap: "wrap", gap: "0 24px" }}>
            <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", width: "100%", marginBottom: 6, fontWeight: 600 }}>Experimental Conditions</div>
            {[
              ["Method", info.method],
              ["Resolution", info.resolution ? `${info.resolution.toFixed(2)} Å` : null],
              ["Crystal pH", info.crystPh],
              ["Crystal Temp", info.crystTemp ? `${info.crystTemp} K` : null],
              ["Collection Temp", info.collectionTemp ? `${info.collectionTemp} K` : null],
              ["R-free", info.rFree ? info.rFree.toFixed(4) : null],
              ["R-work", info.rWork ? info.rWork.toFixed(4) : null],
              ["Organism", info.organism],
              ["Expression", info.expressionSystem],
              ["Deposited", info.deposited ? info.deposited.split("T")[0] : null],
              ["Released", info.released ? info.released.split("T")[0] : null],
            ].filter(([, v]) => v != null).map(([label, value]) => (
              <div key={label} style={{ fontSize: 12, marginBottom: 4, minWidth: 140 }}>
                <span style={{ color: "#556677" }}>{label}: </span>
                <span style={{ color: "#c0c8d8", fontFamily: "'JetBrains Mono', monospace", fontSize: 11 }}>{value}</span>
              </div>
            ))}
            {info.crystDetails && (
              <div style={{ fontSize: 11, color: "#445566", marginTop: 4, width: "100%", lineHeight: 1.4 }}>
                <span style={{ color: "#556677" }}>Growth details: </span>{info.crystDetails}
              </div>
            )}
          </div>
        )}

        {/* ── Functional pH / Enzyme Annotation ── */}
        {info && (info.ecNumber || brendaData || subcellularPH) && (
          <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 16px", marginBottom: 14 }}>
            <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 8, fontWeight: 600 }}>Functional pH Resolution</div>
            <div style={{ display: "flex", flexWrap: "wrap", gap: "8px 24px", alignItems: "flex-start" }}>
              {/* EC Number */}
              {info.ecNumber && (
                <div style={{ fontSize: 12, marginBottom: 2 }}>
                  <span style={{ color: "#556677" }}>EC Number: </span>
                  <span style={{ color: "#4488ff", fontFamily: "'JetBrains Mono', monospace", fontSize: 12, fontWeight: 600 }}>{info.ecNumber}</span>
                </div>
              )}
              {/* Tier 1: BRENDA data */}
              {brendaData && (
                <>
                  <div style={{ fontSize: 12 }}>
                    <span style={{ color: "#556677" }}>pH Optimum: </span>
                    <span style={{ color: "#00cc88", fontFamily: "'JetBrains Mono', monospace", fontSize: 14, fontWeight: 700 }}>{brendaData.ph_median.toFixed(1)}</span>
                    <span style={{ color: "#445566", fontSize: 11 }}> (median, n={brendaData.n_entries})</span>
                  </div>
                  <div style={{ fontSize: 12 }}>
                    <span style={{ color: "#556677" }}>pH Range: </span>
                    <span style={{ color: "#c0c8d8", fontFamily: "'JetBrains Mono', monospace", fontSize: 12 }}>{brendaData.ph_min.toFixed(1)} – {brendaData.ph_max.toFixed(1)}</span>
                  </div>
                  <div style={{ fontSize: 12 }}>
                    <span style={{ color: "#556677" }}>pH Mean: </span>
                    <span style={{ color: "#c0c8d8", fontFamily: "'JetBrains Mono', monospace", fontSize: 12 }}>{brendaData.ph_mean.toFixed(2)}</span>
                  </div>
                  {info.crystPh != null && (
                    <div style={{ fontSize: 12 }}>
                      <span style={{ color: "#556677" }}>ΔpH (cryst − opt): </span>
                      <span style={{
                        color: Math.abs(info.crystPh - brendaData.ph_median) > 1.5 ? "#ff4466" : Math.abs(info.crystPh - brendaData.ph_median) > 0.5 ? "#ffd43b" : "#00cc88",
                        fontFamily: "'JetBrains Mono', monospace", fontSize: 12, fontWeight: 600
                      }}>
                        {(info.crystPh - brendaData.ph_median) >= 0 ? "+" : ""}{(info.crystPh - brendaData.ph_median).toFixed(1)}
                      </span>
                    </div>
                  )}
                  <div style={{ fontSize: 10, padding: "2px 8px", borderRadius: 4, background: "#00cc8815", border: "1px solid #00cc8830", color: "#00cc88", alignSelf: "center" }}>
                    BRENDA SOAP API · {brendaData.organisms?.length || 0} organisms
                  </div>
                </>
              )}
              {/* Tier 2: Subcellular localization pH */}
              {!brendaData && subcellularPH && (
                <>
                  <div style={{ fontSize: 12 }}>
                    <span style={{ color: "#556677" }}>Functional pH: </span>
                    <span style={{ color: "#da77f2", fontFamily: "'JetBrains Mono', monospace", fontSize: 14, fontWeight: 700 }}>{subcellularPH.ph.toFixed(1)}</span>
                  </div>
                  <div style={{ fontSize: 12 }}>
                    <span style={{ color: "#556677" }}>Compartment: </span>
                    <span style={{ color: "#c0c8d8", fontSize: 12 }}>{subcellularPH.label}</span>
                  </div>
                  <div style={{ fontSize: 12 }}>
                    <span style={{ color: "#556677" }}>pH Range: </span>
                    <span style={{ color: "#c0c8d8", fontFamily: "'JetBrains Mono', monospace", fontSize: 12 }}>{subcellularPH.range[0].toFixed(1)} – {subcellularPH.range[1].toFixed(1)}</span>
                  </div>
                  {info.crystPh != null && (
                    <div style={{ fontSize: 12 }}>
                      <span style={{ color: "#556677" }}>ΔpH (cryst − func): </span>
                      <span style={{
                        color: Math.abs(info.crystPh - subcellularPH.ph) > 1.5 ? "#ff4466" : Math.abs(info.crystPh - subcellularPH.ph) > 0.5 ? "#ffd43b" : "#00cc88",
                        fontFamily: "'JetBrains Mono', monospace", fontSize: 12, fontWeight: 600
                      }}>
                        {(info.crystPh - subcellularPH.ph) >= 0 ? "+" : ""}{(info.crystPh - subcellularPH.ph).toFixed(1)}
                      </span>
                    </div>
                  )}
                  <div style={{ fontSize: 10, padding: "2px 8px", borderRadius: 4, background: "#da77f215", border: "1px solid #da77f230", color: "#da77f2", alignSelf: "center" }}>
                    UniProt subcellular localization
                  </div>
                </>
              )}
              {/* No data yet */}
              {info.ecNumber && !brendaData && !subcellularPH && (
                <div style={{ fontSize: 11, color: "#445566", fontStyle: "italic" }}>
                  Looking up functional pH...
                </div>
              )}
            </div>
            {/* Organisms list for BRENDA */}
            {brendaData?.organisms?.length > 0 && (
              <details style={{ marginTop: 6 }}>
                <summary style={{ fontSize: 10, color: "#556677", cursor: "pointer" }}>Organisms ({brendaData.organisms.length})</summary>
                <div style={{ fontSize: 11, color: "#445566", marginTop: 4, lineHeight: 1.6, maxHeight: 80, overflowY: "auto" }}>
                  {brendaData.organisms.join(" · ")}
                </div>
              </details>
            )}
            {/* All locations for subcellular */}
            {!brendaData && subcellularPH?.allLocations?.length > 1 && (
              <details style={{ marginTop: 6 }}>
                <summary style={{ fontSize: 10, color: "#556677", cursor: "pointer" }}>All locations ({subcellularPH.allLocations.length})</summary>
                <div style={{ fontSize: 11, color: "#445566", marginTop: 4, lineHeight: 1.6 }}>
                  {subcellularPH.allLocations.join(" · ")}
                </div>
              </details>
            )}
          </div>
        )}

        {/* ── BRENDA Credentials Banner ── */}
        {brendaNeedsCredentials && info?.ecNumber && !brendaData && (
          <BrendaCredentialsBanner onOpenSettings={() => setShowSettings(true)} />
        )}

        {/* ── Metrics ── */}
        <div style={{ display: "grid", gridTemplateColumns: "repeat(6, 1fr)", gap: 10, marginBottom: 18 }}>
          {[
            ["PDB ID", pdbId || filePdbId || "—", "#4488ff"],
            ["Residues", info?.sequence?.length || "—", "#d0d8e8"],
            ["pI", pI ? pI.toFixed(2) : "—", "#00cc88"],
            ["Net Charge", `${charge >= 0 ? "+" : ""}${charge.toFixed(1)}`, charge > 0 ? "#4dabf7" : "#ff6b6b"],
            ["pH-Sensitive", nSensitive, "#ff4466"],
            ["Titratable", titratable.length, "#ffd43b"],
          ].map(([label, value, color]) => (
            <div key={label} style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "10px 14px" }}>
              <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.05em", marginBottom: 2 }}>{label}</div>
              <div style={{ fontSize: 18, fontWeight: 700, color, fontFamily: "'JetBrains Mono', monospace" }}>{value}</div>
            </div>
          ))}
        </div>

        {/* ── pKa Data Source Badge ── */}
        {pkaSource && (
          <div style={{ marginBottom: 12, display: "flex", alignItems: "center", gap: 8 }}>
            <span style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.05em" }}>pKa Source:</span>
            <span style={{
              fontSize: 11, fontWeight: 600, padding: "2px 8px", borderRadius: 4,
              background: pkaSource === "propka" ? "#00cc8820" : pkaSource === "standard_fallback" ? "#ffd43b20" : "#ff446620",
              color: pkaSource === "propka" ? "#00cc88" : pkaSource === "standard_fallback" ? "#ffd43b" : "#ff4466",
              border: `1px solid ${pkaSource === "propka" ? "#00cc8840" : pkaSource === "standard_fallback" ? "#ffd43b40" : "#ff446640"}`,
            }}>
              {pkaSource === "propka" ? "PROPKA3 (structure-based)" : pkaSource === "standard_fallback" ? "Standard pKa (no PROPKA)" : "Approximation (random perturbation)"}
            </span>
          </div>
        )}

        {/* ── pH Slider ── */}
        <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "14px 18px", marginBottom: 18, display: "flex", alignItems: "center", gap: 16 }}>
          <div style={{ fontSize: 11, color: "#556677", textTransform: "uppercase", letterSpacing: "0.05em", whiteSpace: "nowrap" }}>Target pH</div>
          <input type="range" min="0" max="14" step="0.1" value={ph} onChange={e => setPh(parseFloat(e.target.value))} style={{ flex: 1 }} />
          <div style={{ fontSize: 22, fontWeight: 700, fontFamily: "'JetBrains Mono', monospace", color: "#ff4466", minWidth: 50, textAlign: "right" }}>{ph.toFixed(1)}</div>
        </div>

        {/* ── 3D Viewer Controls + Viewer ── */}
        <div style={{ display: "grid", gridTemplateColumns: "220px 1fr", gap: 14, marginBottom: 16 }}>

          {/* ── Left Controls Panel ── */}
          <div style={{ display: "flex", flexDirection: "column", gap: 10 }}>

            {/* Visualization Layers */}
            <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 14px" }}>
              <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 8, fontWeight: 600 }}>Visualization Layers</div>
              {[
                ["ph_sensitive", "pH-Sensitive (±1.0)", "#ff4466"],
                ["titratable", "Titratable Residues", "#00cc88"],
                ["protonation", "Protonation States", "#4dabf7"],
                ["pka_heatmap", "pKa Heatmap", "#ffd43b"],
              ].map(([key, label, color]) => (
                <label key={key} style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer", padding: "4px 0", fontSize: 12 }}>
                  <div onClick={() => toggleLayer(key)} style={{
                    width: 16, height: 16, borderRadius: 4, border: `2px solid ${layers[key] ? color : "#334"}`,
                    background: layers[key] ? `${color}30` : "transparent", display: "flex", alignItems: "center", justifyContent: "center",
                    transition: "all 0.15s", cursor: "pointer",
                  }}>
                    {layers[key] && <span style={{ color, fontSize: 11, lineHeight: 1 }}>✓</span>}
                  </div>
                  <span style={{ color: layers[key] ? "#d0d8e8" : "#556677" }}>{label}</span>
                </label>
              ))}
            </div>

            {/* Sequence Range */}
            <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 14px" }}>
              <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 8, fontWeight: 600 }}>Sequence Range</div>
              <div style={{ display: "flex", gap: 6, marginBottom: 8 }}>
                <input type="number" value={rangeStart} onChange={e => setRangeStart(e.target.value)} placeholder="From" min="1" style={{ width: "50%", padding: "6px 8px", background: "#080c14", border: "1px solid #1e2a3a", borderRadius: 4, color: "#d0d8e8", fontSize: 12, fontFamily: "'JetBrains Mono', monospace" }} />
                <input type="number" value={rangeEnd} onChange={e => setRangeEnd(e.target.value)} placeholder="To" min="1" style={{ width: "50%", padding: "6px 8px", background: "#080c14", border: "1px solid #1e2a3a", borderRadius: 4, color: "#d0d8e8", fontSize: 12, fontFamily: "'JetBrains Mono', monospace" }} />
              </div>
              <div style={{ display: "flex", gap: 6 }}>
                <button onClick={applyRangeFocus} style={{ flex: 1, padding: "6px 0", background: "#00cc8820", border: "1px solid #00cc8840", borderRadius: 4, color: "#00cc88", fontSize: 11, fontWeight: 600, cursor: "pointer" }}>Focus</button>
                <button onClick={clearRange} style={{ flex: 1, padding: "6px 0", background: "#33445520", border: "1px solid #33445540", borderRadius: 4, color: "#8899aa", fontSize: 11, fontWeight: 600, cursor: "pointer" }}>Clear</button>
              </div>
              {info?.sequence && (
                <div style={{ fontSize: 10, color: "#334455", marginTop: 6 }}>Range: 1 – {info.sequence.length}</div>
              )}
            </div>

            {/* Quick Select */}
            <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "12px 14px" }}>
              <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 8, fontWeight: 600 }}>Quick Select</div>
              <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
                <button onClick={quickSelectTitratable} style={{ padding: "7px 0", background: "#0d1520", border: "1px solid #00cc8840", borderRadius: 4, color: "#00cc88", fontSize: 11, fontWeight: 600, cursor: "pointer" }}>Titratable ({titratable.length})</button>
                <button onClick={quickSelectPhSensitive} style={{ padding: "7px 0", background: "#0d1520", border: "1px solid #ff446640", borderRadius: 4, color: "#ff4466", fontSize: 11, fontWeight: 600, cursor: "pointer" }}>pH-Sensitive ({titratable.filter(r => Math.abs(r.pka - ph) < 1).length})</button>
                <button onClick={quickClearSelection} style={{ padding: "7px 0", background: "#0d1520", border: "1px solid #33445540", borderRadius: 4, color: "#8899aa", fontSize: 11, fontWeight: 600, cursor: "pointer" }}>Clear Selection</button>
              </div>
            </div>

            {/* Legend */}
            <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: 8, padding: "10px 14px" }}>
              <div style={{ fontSize: 10, color: "#556677", textTransform: "uppercase", letterSpacing: "0.06em", marginBottom: 6, fontWeight: 600 }}>Legend</div>
              {[
                ["#ff4466", "pH-Sensitive (±1.0)"],
                ["#ff6b6b", "Acidic (ASP, GLU)"],
                ["#4dabf7", "Basic (HIS, LYS)"],
                ["#ffd43b", "Thiol (CYS)"],
                ["#da77f2", "Phenolic (TYR)"],
                ["#4dabf7", "Protonated"],
                ["#ff6b6b", "Deprotonated"],
              ].map(([color, label]) => (
                <div key={label} style={{ display: "flex", alignItems: "center", gap: 6, padding: "2px 0", fontSize: 11 }}>
                  <div style={{ width: 10, height: 10, borderRadius: 2, background: color, flexShrink: 0 }} />
                  <span style={{ color: "#8899aa" }}>{label}</span>
                </div>
              ))}
            </div>
          </div>

          {/* ── Right: 3D Viewer ── */}
          <div>
            <MolstarViewer ref={molstarRef} pdbId={pdbId} customUrl={customUrl} customData={customData} customFormat={customFormat} />
            <div style={{ fontSize: 11, color: "#667788", marginTop: 6 }}>
              Mol* Viewer — hover to highlight · click to focus · sequence panel above 3D view · scroll to zoom
            </div>
          </div>
        </div>

        {/* ── Analysis Tabs ── */}
        <div style={{ display: "flex", gap: 2, marginBottom: 2 }}>
          {[["titration", "Titration Curve"], ["pka", "pKa Along Sequence"], ["table", "Titratable Residues"], ["sequence", "Full Sequence"]].map(([key, label]) => (
            <button key={key} onClick={() => setActiveTab(key)} style={{ padding: "10px 20px", background: activeTab === key ? "#0d1520" : "transparent", border: `1px solid ${activeTab === key ? "#1e2a3a" : "transparent"}`, borderBottom: activeTab === key ? "1px solid #0d1520" : "1px solid #1e2a3a", borderRadius: "8px 8px 0 0", color: activeTab === key ? "#d0d8e8" : "#556677", fontSize: 13, fontWeight: activeTab === key ? 600 : 400, cursor: "pointer" }}>{label}</button>
          ))}
        </div>
        <div style={{ background: "#0d1520", border: "1px solid #1e2a3a", borderRadius: "0 8px 8px 8px", padding: 20, marginBottom: 24 }}>
          {activeTab === "titration" && (
            <div>
              <div style={{ fontSize: 11, color: "#556677", marginBottom: 10 }}>Click anywhere on the curve to set pH · pI and crystallization pH marked</div>
              <TitrationCurve curve={curve} ph={ph} pI={pI} crystPh={info?.crystPh} brendaPh={brendaData?.ph_median ?? subcellularPH?.ph} onPhChange={setPh} />
            </div>
          )}
          {activeTab === "pka" && (
            <div>
              <div style={{ fontSize: 11, color: "#556677", marginBottom: 10 }}>Bars colored by proximity to current pH: <span style={{color:"#ff4466"}}>●</span> within ±1 <span style={{color:"#ffcc00"}}>●</span> within ±2 <span style={{color:"#4488ff"}}>●</span> far</div>
              <PkaBarChart titratable={titratable} ph={ph} />
            </div>
          )}
          {activeTab === "table" && (
            <div>
              <div style={{ display: "flex", gap: 10, marginBottom: 12 }}>
                {[["Acidic (ASP, GLU)", nAcidic, "#ff6b6b"], ["Basic (HIS, LYS)", nBasic, "#4dabf7"], ["Thiol (CYS)", titratable.filter(r => r.cls === "Thiol").length, "#ffd43b"], ["Phenolic (TYR)", titratable.filter(r => r.cls === "Phenolic").length, "#da77f2"]].map(([label, n, color]) => (
                  <div key={label} style={{ background: `${color}10`, border: `1px solid ${color}30`, borderRadius: 6, padding: "6px 14px", fontSize: 12 }}>
                    <span style={{ color }}>{n}</span> <span style={{ color: "#8899aa" }}>{label}</span>
                  </div>
                ))}
              </div>
              <TitratableTable titratable={titratable.map(r => ({ ...r, delta: r.pka - ph }))} ph={ph} pkaSource={pkaSource} />
            </div>
          )}
          {activeTab === "sequence" && info?.sequence && (
            <div>
              <div style={{ fontSize: 11, color: "#556677", marginBottom: 10 }}>
                {info.sequence.length} residues · Titratable residues highlighted · pKa-sensitive residues (±1.0 of current pH) shown in red
              </div>
              <FullSequenceTable sequence={info.sequence} titratable={titratable} ph={ph} />
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
