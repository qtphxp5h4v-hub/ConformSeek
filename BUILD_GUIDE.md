# ConformSeek — Desktop App Build Guide

pH-Dependent Protein Conformational Analysis

## Architecture

```
ConformSeek.app / ConformSeek.exe
├── Electron shell (native window)
├── React frontend (build/ → static HTML/JS/CSS)
└── Python backend (PyInstaller binary → conformseek_api)
    ├── FastAPI + Uvicorn (HTTP server on localhost:8000)
    ├── PROPKA3 (pKa prediction)
    ├── BioPython (structure parsing)
    └── Zeep (BRENDA SOAP API client)
```

When the app launches:
1. Electron spawns the Python backend binary
2. A splash screen shows while waiting for the backend health check
3. Once healthy, the React UI loads and communicates with the backend over localhost:8000

---

## Prerequisites

- **Node.js** ≥ 18 (for React + Electron)
- **Python** ≥ 3.9 (for PyInstaller compilation)
- **pip** packages: `pyinstaller fastapi uvicorn propka biopython numpy requests zeep lxml`

---

## Project Structure

```
conformseek-viewer/
├── main.js                    ← Electron main process
├── package.json               ← npm scripts + electron-builder config
├── assets/
│   ├── entitlements.mac.plist ← macOS entitlements for code signing
│   ├── icon.icns              ← macOS app icon (you create this)
│   └── icon.ico               ← Windows app icon (you create this)
├── public/
│   └── index.html             ← React HTML shell (from create-react-app)
├── src/
│   ├── ConformSeekApp.jsx     ← Main React component
│   ├── index.js               ← React entry point
│   └── index.css              ← Base styles
├── conformseek_server.py      ← Python backend (source, for PyInstaller)
├── python-dist/               ← PyInstaller output goes here
│   └── conformseek_api        ← Compiled Python binary
├── build/                     ← React production build (generated)
└── dist/                      ← Final .dmg / .exe (generated)
```

---

## Step-by-Step Build

### 1. Set Up the React Project (one-time)

If you're starting from your existing `conformseek-viewer` folder:

```bash
cd /path/to/conformseek-viewer

# Copy the new Electron files into your project root
cp /path/to/main.js .
cp /path/to/package.json .     # ← merge with your existing one (see note)
mkdir -p assets
cp /path/to/assets/entitlements.mac.plist assets/

# Install Electron dev dependencies
npm install --save-dev electron electron-builder concurrently wait-on
```

> **Note on package.json**: If you already have a `package.json` from `create-react-app`,
> you need to merge in the new fields. The key additions are:
> - `"main": "main.js"` (top level)
> - `"homepage": "./"` (top level — critical for Electron file:// loading)
> - The `electron:*` scripts
> - The `"build"` section (electron-builder config)
> - The devDependencies (electron, electron-builder, concurrently, wait-on)

### 2. Compile the Python Backend with PyInstaller

```bash
# Make sure you're in a Python environment with all dependencies
pip install pyinstaller fastapi uvicorn propka biopython numpy requests zeep lxml

# Run PyInstaller (from your project root where conformseek_server.py lives)
pyinstaller \
  --onefile \
  --name conformseek_api \
  --collect-all propka \
  --collect-all Bio \
  --collect-all zeep \
  --collect-all lxml \
  conformseek_server.py

# Move the binary to where Electron expects it
mkdir -p python-dist
cp dist/conformseek_api python-dist/
```

**Important**: PyInstaller compiles for the platform you run it on.
- On your Mac (M4 Max) → produces an ARM64 macOS binary
- On Windows → produces a .exe
- You cannot cross-compile. Use GitHub Actions for multi-platform builds (see below).

### 3. Test in Development Mode

This runs the React dev server (port 3000) and Electron side by side:

```bash
# Terminal 1: Start the Python backend directly (for faster iteration)
python conformseek_server.py

# Terminal 2: Start Electron + React dev server
npm run electron:dev
```

Or on Windows:
```bash
npm run electron:dev:win
```

You should see:
1. A splash screen saying "Starting analysis server..."
2. The ConformSeek UI in a native window
3. Full functionality (PROPKA, BRENDA, Mol* viewer, etc.)

### 4. Build the Production App

```bash
# Build React → build/ folder
# Then package with electron-builder → dist/ folder
npm run electron:build:mac    # → dist/ConformSeek-2.0.0-arm64.dmg
npm run electron:build:win    # → dist/ConformSeek Setup 2.0.0.exe
npm run electron:build:linux  # → dist/ConformSeek-2.0.0.AppImage
```

### 5. Run the .dmg

1. Open `dist/ConformSeek-2.0.0-arm64.dmg`
2. Drag ConformSeek into Applications
3. **First launch**: Right-click → Open (macOS Gatekeeper will warn about unsigned app)
4. After that, it opens normally from Launchpad/Spotlight

---

## App Icons (optional but recommended)

Create a 1024x1024 PNG icon, then convert:

**macOS (.icns)**:
```bash
mkdir icon.iconset
sips -z 16 16     icon.png --out icon.iconset/icon_16x16.png
sips -z 32 32     icon.png --out icon.iconset/icon_16x16@2x.png
sips -z 32 32     icon.png --out icon.iconset/icon_32x32.png
sips -z 64 64     icon.png --out icon.iconset/icon_32x32@2x.png
sips -z 128 128   icon.png --out icon.iconset/icon_128x128.png
sips -z 256 256   icon.png --out icon.iconset/icon_128x128@2x.png
sips -z 256 256   icon.png --out icon.iconset/icon_256x256.png
sips -z 512 512   icon.png --out icon.iconset/icon_256x256@2x.png
sips -z 512 512   icon.png --out icon.iconset/icon_512x512.png
sips -z 1024 1024 icon.png --out icon.iconset/icon_512x512@2x.png
iconutil -c icns icon.iconset -o assets/icon.icns
```

**Windows (.ico)**: Use an online converter or ImageMagick:
```bash
convert icon.png -define icon:auto-resize=256,128,64,48,32,16 assets/icon.ico
```

If you don't create icons, electron-builder will use a default Electron icon — the app still works fine.

---

## GitHub Distribution

### Release Structure

```
ConformSeek/
├── README.md
├── LICENSE
├── releases/
│   ├── ConformSeek-2.0.0-arm64.dmg      (macOS Apple Silicon)
│   ├── ConformSeek-2.0.0-x64.dmg        (macOS Intel)
│   ├── ConformSeek-Setup-2.0.0.exe       (Windows x64)
│   └── ConformSeek-2.0.0.AppImage        (Linux x64)
└── src/                                   (source code)
```

### GitHub Actions (automated cross-platform builds)

Create `.github/workflows/build.yml`:

```yaml
name: Build ConformSeek

on:
  push:
    tags: ['v*']

jobs:
  build-python:
    strategy:
      matrix:
        os: [macos-14, windows-latest, ubuntu-latest]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.11'
      - run: pip install pyinstaller fastapi uvicorn propka biopython numpy requests zeep lxml
      - run: pyinstaller --onefile --name conformseek_api --collect-all propka --collect-all Bio --collect-all zeep --collect-all lxml conformseek_server.py
      - uses: actions/upload-artifact@v4
        with:
          name: python-binary-${{ matrix.os }}
          path: dist/conformseek_api*

  build-electron:
    needs: build-python
    strategy:
      matrix:
        include:
          - os: macos-14
            electron_target: mac
          - os: windows-latest
            electron_target: win
          - os: ubuntu-latest
            electron_target: linux
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with:
          node-version: '20'
      - uses: actions/download-artifact@v4
        with:
          name: python-binary-${{ matrix.os }}
          path: python-dist/
      - run: chmod +x python-dist/conformseek_api || true
      - run: npm ci
      - run: npm run build
      - run: npx electron-builder --${{ matrix.electron_target }}
      - uses: actions/upload-artifact@v4
        with:
          name: release-${{ matrix.os }}
          path: |
            dist/*.dmg
            dist/*.exe
            dist/*.AppImage
```

Push a tag to trigger: `git tag v2.0.0 && git push --tags`

---

## Troubleshooting

**"Python backend not starting"**
- Check that `python-dist/conformseek_api` exists and is executable (`chmod +x`)
- Test it manually: `./python-dist/conformseek_api` — should start a server on port 8000
- On macOS, you may need to allow it in System Preferences → Privacy & Security

**"App shows blank window"**
- In dev: make sure `npm start` is running (React dev server on port 3000)
- In production: make sure `npm run build` ran before `electron-builder`

**"BRENDA queries fail"**
- Open Settings (⚙) and check credential status
- Make sure zeep was included in the PyInstaller build (`--collect-all zeep`)

**Binary size is large (200-400MB)**
- This is normal — PyInstaller bundles the entire Python runtime + NumPy + BioPython
- UPX compression can reduce it: `pyinstaller --onefile --upx-dir=/path/to/upx ...`

**macOS Gatekeeper blocks the app**
- Without code signing: Right-click → Open → Open (first launch only)
- With signing: Obtain an Apple Developer certificate and add to electron-builder config
