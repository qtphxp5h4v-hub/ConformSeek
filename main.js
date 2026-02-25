// main.js — Electron main process for ConformSeek
// Launches the Python backend server, waits for it to be healthy,
// then shows the React UI in a native window.

const { app, BrowserWindow, dialog } = require("electron");
const path = require("path");
const fs = require("fs");
const { spawn, execSync } = require("child_process");
const http = require("http");

let mainWindow;
let pythonProcess;
let splashWindow;

const PYTHON_PORT = 8000;
const HEALTH_URL = `http://localhost:${PYTHON_PORT}/health`;

// ── Resolve path to the PyInstaller binary ──
function getPythonExecutablePath() {
  const isDev = !app.isPackaged;
  const ext = process.platform === "win32" ? ".exe" : "";
  const binaryName = `conformseek_api${ext}`;

  if (isDev) {
    return path.join(__dirname, "python-dist", binaryName);
  } else {
    return path.join(process.resourcesPath, binaryName);
  }
}

// ── Kill any stale conformseek_api processes and free port 8000 ──
function cleanupStaleProcesses() {
  try {
    if (process.platform === "win32") {
      // Windows: kill by name
      execSync("taskkill /F /IM conformseek_api.exe 2>nul", { stdio: "ignore" });
    } else {
      // macOS/Linux: kill any existing conformseek_api processes
      execSync("pkill -9 -f conformseek_api 2>/dev/null || true", { stdio: "ignore" });

      // Also kill anything holding port 8000
      try {
        const lsofOutput = execSync(`lsof -ti :${PYTHON_PORT} 2>/dev/null`).toString().trim();
        if (lsofOutput) {
          const pids = lsofOutput.split("\n").filter(Boolean);
          for (const pid of pids) {
            console.log(`[ConformSeek] Killing stale process on port ${PYTHON_PORT}: PID ${pid}`);
            execSync(`kill -9 ${pid} 2>/dev/null || true`, { stdio: "ignore" });
          }
        }
      } catch (e) {
        // No process on port — that's fine
      }
    }
    // Brief pause to let the OS release the port
    return new Promise((resolve) => setTimeout(resolve, 500));
  } catch (e) {
    console.log("[ConformSeek] Cleanup (non-critical):", e.message);
    return Promise.resolve();
  }
}

// ── Start the Python backend ──
function startPythonBackend() {
  const exePath = getPythonExecutablePath();
  console.log(`[ConformSeek] Binary path: ${exePath}`);

  // Verify the binary exists
  if (!fs.existsSync(exePath)) {
    const msg = `Python backend binary not found at:\n${exePath}`;
    console.error(`[ConformSeek] ${msg}`);
    dialog.showErrorBox("ConformSeek — Missing Backend", msg);
    return false;
  }

  // Verify it's executable (macOS/Linux)
  if (process.platform !== "win32") {
    try {
      fs.accessSync(exePath, fs.constants.X_OK);
    } catch (e) {
      console.log("[ConformSeek] Setting executable permission on binary...");
      try {
        fs.chmodSync(exePath, 0o755);
      } catch (chmodErr) {
        console.error("[ConformSeek] Failed to set executable permission:", chmodErr.message);
      }
    }
  }

  console.log(`[ConformSeek] Spawning Python backend...`);

  pythonProcess = spawn(exePath, [], {
    stdio: ["ignore", "pipe", "pipe"],
    env: { ...process.env, PYTHONUNBUFFERED: "1" },
    detached: false,
  });

  pythonProcess.stdout.on("data", (data) => {
    const lines = data.toString().trim();
    console.log(`[Python] ${lines}`);
  });

  pythonProcess.stderr.on("data", (data) => {
    const lines = data.toString().trim();
    // Uvicorn logs to stderr by default — not necessarily errors
    console.log(`[Python] ${lines}`);
  });

  pythonProcess.on("error", (err) => {
    console.error("[ConformSeek] Failed to spawn Python backend:", err.message);
    dialog.showErrorBox(
      "ConformSeek — Backend Error",
      `Could not start the Python backend server.\n\n${err.message}\n\nPath: ${exePath}\n\nIf macOS is blocking it, go to:\nSystem Settings → Privacy & Security → Allow Anyway`
    );
    pythonProcess = null;
    return false;
  });

  pythonProcess.on("exit", (code, signal) => {
    console.log(`[ConformSeek] Python backend exited (code=${code}, signal=${signal})`);
    // If it crashed immediately, it might be Gatekeeper blocking it
    if (code !== null && code !== 0 && mainWindow) {
      console.error(`[ConformSeek] Backend crashed with exit code ${code}`);
    }
    pythonProcess = null;
  });

  console.log(`[ConformSeek] Python process spawned (PID: ${pythonProcess.pid})`);
  return true;
}

// ── Wait for the backend to become healthy ──
function waitForBackend(maxRetries = 40, intervalMs = 500) {
  return new Promise((resolve, reject) => {
    let attempts = 0;

    const check = () => {
      attempts++;
      if (attempts % 5 === 0) {
        console.log(`[ConformSeek] Health check attempt ${attempts}/${maxRetries}...`);
      }

      // If the process already died, don't keep trying
      if (!pythonProcess) {
        reject(new Error("Python backend process is not running"));
        return;
      }

      const req = http.get(HEALTH_URL, (res) => {
        let body = "";
        res.on("data", (chunk) => { body += chunk; });
        res.on("end", () => {
          if (res.statusCode === 200) {
            console.log(`[ConformSeek] Backend healthy after ${attempts} attempts`);
            try {
              const health = JSON.parse(body);
              console.log(`[ConformSeek] Backend status:`, JSON.stringify(health));
            } catch (e) { /* ignore parse errors */ }
            resolve();
          } else {
            retry();
          }
        });
      });

      req.on("error", () => retry());
      req.setTimeout(2000, () => {
        req.destroy();
        retry();
      });
    };

    const retry = () => {
      if (attempts >= maxRetries) {
        reject(new Error(`Backend did not become healthy after ${maxRetries} attempts (${maxRetries * intervalMs / 1000}s)`));
      } else {
        setTimeout(check, intervalMs);
      }
    };

    check();
  });
}

// ── Force kill the Python backend ──
function killPythonBackend() {
  if (pythonProcess) {
    const pid = pythonProcess.pid;
    console.log(`[ConformSeek] Stopping Python backend (PID: ${pid})...`);

    try {
      // Try graceful shutdown first
      pythonProcess.kill("SIGTERM");
    } catch (e) { /* ignore */ }

    // Force kill after 2 seconds if still running
    setTimeout(() => {
      try {
        process.kill(pid, 0); // Check if still alive (throws if dead)
        console.log(`[ConformSeek] Force killing Python backend (PID: ${pid})...`);
        process.kill(pid, "SIGKILL");
      } catch (e) {
        // Process already dead — good
      }
    }, 2000);

    pythonProcess = null;
  }

  // Also clean up any orphaned processes on the port
  if (process.platform !== "win32") {
    try {
      execSync(`lsof -ti :${PYTHON_PORT} | xargs kill -9 2>/dev/null || true`, { stdio: "ignore" });
    } catch (e) { /* ignore */ }
  }
}

// ── Create the splash/loading window ──
function createSplashWindow() {
  splashWindow = new BrowserWindow({
    width: 340,
    height: 220,
    frame: false,
    transparent: true,
    resizable: false,
    alwaysOnTop: true,
    webPreferences: { nodeIntegration: false, contextIsolation: true },
  });

  splashWindow.loadURL(
    `data:text/html;charset=utf-8,${encodeURIComponent(`
    <!DOCTYPE html>
    <html>
    <head>
      <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          background: #0d1520;
          color: #d0d8e8;
          display: flex;
          align-items: center;
          justify-content: center;
          height: 100vh;
          border-radius: 12px;
          border: 1px solid #1e2a3a;
          overflow: hidden;
        }
        .container { text-align: center; }
        .logo {
          width: 48px; height: 48px; border-radius: 12px;
          background: linear-gradient(135deg, #4488ff, #00cc88);
          display: inline-flex; align-items: center; justify-content: center;
          font-size: 24px; margin-bottom: 14px;
        }
        h1 { font-size: 18px; font-weight: 700; margin-bottom: 6px; letter-spacing: -0.3px; }
        p { font-size: 12px; color: #556677; margin-bottom: 16px; }
        .spinner {
          width: 24px; height: 24px;
          border: 3px solid #1e2a3a; border-top-color: #4488ff;
          border-radius: 50%; margin: 0 auto;
          animation: spin 1s linear infinite;
        }
        @keyframes spin { to { transform: rotate(360deg); } }
      </style>
    </head>
    <body>
      <div class="container">
        <div class="logo">🔬</div>
        <h1>ConformSeek</h1>
        <p>Starting analysis server...</p>
        <div class="spinner"></div>
      </div>
    </body>
    </html>
  `)}`
  );
}

// ── Create the main application window ──
function createMainWindow() {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 960,
    minHeight: 600,
    show: false,
    backgroundColor: "#080c14",
    titleBarStyle: "hiddenInset",
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
    },
  });

  const isDev = !app.isPackaged;
  if (isDev) {
    mainWindow.loadURL("http://localhost:3000");
  } else {
    mainWindow.loadFile(path.join(__dirname, "build", "index.html"));
  }

  mainWindow.once("ready-to-show", () => {
    if (splashWindow && !splashWindow.isDestroyed()) {
      splashWindow.close();
      splashWindow = null;
    }
    mainWindow.show();
  });

  mainWindow.on("closed", () => {
    mainWindow = null;
  });
}

// ── App lifecycle ──
app.whenReady().then(async () => {
  console.log("[ConformSeek] App starting...");
  console.log(`[ConformSeek] Packaged: ${app.isPackaged}`);
  console.log(`[ConformSeek] Platform: ${process.platform} (${process.arch})`);
  console.log(`[ConformSeek] Resources: ${process.resourcesPath}`);

  createSplashWindow();

  // Kill any stale processes from previous runs
  await cleanupStaleProcesses();

  // Start the Python backend
  const started = startPythonBackend();

  if (started) {
    try {
      await waitForBackend(40, 500); // up to 20 seconds
    } catch (err) {
      console.error("[ConformSeek]", err.message);
      dialog.showMessageBox({
        type: "warning",
        title: "ConformSeek — Backend Slow",
        message: "The analysis server is taking longer than expected to start.",
        detail: "The app will open, but PROPKA and BRENDA features may not work until the server finishes loading.\n\nIf this persists, check System Settings → Privacy & Security to make sure conformseek_api is allowed.",
        buttons: ["OK"],
      });
    }
  }

  createMainWindow();
});

// Kill the Python server when the app quits
app.on("will-quit", () => {
  killPythonBackend();
});

// Also handle before-quit for macOS Cmd+Q
app.on("before-quit", () => {
  killPythonBackend();
});

// macOS: quit when all windows closed
app.on("window-all-closed", () => {
  if (process.platform !== "darwin") {
    app.quit();
  }
});

// macOS: re-create window when dock icon clicked
app.on("activate", () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createMainWindow();
  }
});
