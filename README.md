# ConformSeek
ConformSeek aims to resolve a major issue in modern biochemical and molecular analysis by generating predictions to determine how variable pH environments will affect protein conformation. 

---

## Installation & Setup

Because this app is currently in independent development, your operating system might flag it as "unverified." Don't worry—this is just because the app hasn't been "signed" with a paid developer certificate yet.

### macOS

1. Download the `.dmg` or `.app` file.
2. Move **ConformSeek** to your **Applications** folder.
3. If you see a "Developer cannot be verified" or "App is damaged" error, open your **Terminal** and run:
```bash
xattr -cr /Applications/ConformSeek.app

```


4. Right-click the app and select **Open**.

### Windows

1. Download the `.exe` installer.
2. Run the installer. You will likely see a blue **"Windows protected your PC"** screen.
3. Click **"More info"**.
4. Click **"Run anyway"** at the bottom of the window.

### Linux

1. Download the `.AppImage` file.
2. Right-click the file, go to **Properties > Permissions**, and check **"Allow executing file as program."**
3. Alternatively, run this in your terminal:
```bash
chmod +x ConformSeek.AppImage
./ConformSeek.AppImage

```


*Note: If the app fails to launch on some distributions, try running it with `./ConformSeek.AppImage --no-sandbox`.*

---

## Built With

* [React](https://reactjs.org/) - Frontend UI
* [Python](https://www.python.org/) - Backend server for API requests
* [Electron](https://www.electronjs.org/) - Desktop Framework
* [Electron Builder](https://www.electron.build/) - Packaging and Distribution

## Development

If you want to run this project locally:

1. **Clone the repo:**
```bash
git clone https://github.com/qtphxp5h4v-hub/ConformSeek.git

```


2. **Install dependencies:**
```bash
npm install

```


3. **Start the app:**
```bash
npm start

```



---

### Why the security warnings?

> [!NOTE]
> To remove these warnings, developers must pay annual fees to Apple ($99/yr) and Windows ($200+/yr). To keep this project free and open-source, I have chosen to skip the formal "Signing" process for now. The code is fully open for you to audit!

---
