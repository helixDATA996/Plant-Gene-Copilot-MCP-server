# Agent Human Evaluation System Runtime Package

## Contents

- `index.html`: Evaluation webpage. The default interface language is English, with a Chinese toggle.
- `server.py`: Flask backend for serving the webpage, validating submissions, preventing duplicates, and writing to MySQL.
- `images/`: All image resources used by the webpage.
- `requirements.txt`: Python dependencies.
- `start_evaluation_system.bat`: Windows script that installs dependencies and starts the server.

Historical CSV, JSON, DOCX, and other experiment files are excluded from this runtime package.

## Prerequisites

1. Install Python 3.10 or later.
2. Install and start the MySQL service.
3. Make sure the MySQL user can create databases and tables.

## Start the System

Double-click `start_evaluation_system.bat`, then open:

`http://localhost:5000`

The backend automatically creates the `agent_eval` database and `scores` table during the first health check or submission.

## MySQL Configuration

The default configuration is:

- Host: `localhost`
- Port: `3306`
- User: `root`
- Password: `12345678`
- Database: `agent_eval`

If your configuration is different, set the environment variables in PowerShell before starting:

```powershell
$env:MYSQL_HOST = "localhost"
$env:MYSQL_PORT = "3306"
$env:MYSQL_USER = "root"
$env:MYSQL_PASSWORD = "your-mysql-password"
$env:AGENT_EVAL_DB = "agent_eval"
py server.py
```

Do not commit a real database password to a repository or share it with evaluators.

## LAN Access

For local evaluation, use `http://localhost:5000`.

To allow access from another computer on the same LAN, use the host computer's LAN IP, for example `http://192.168.1.10:5000`. Allow TCP port 5000 through Windows Firewall and expose MySQL only to trusted networks.
