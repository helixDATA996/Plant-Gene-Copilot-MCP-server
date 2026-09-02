from datetime import datetime
import os
import uuid

import pymysql
from flask import Flask, jsonify, request, send_from_directory


app = Flask(__name__)
DIR = os.path.dirname(os.path.abspath(__file__))

DB_NAME = os.getenv("AGENT_EVAL_DB", "agent_eval")
DB_CONFIG = {
    "host": os.getenv("MYSQL_HOST", "localhost"),
    "port": int(os.getenv("MYSQL_PORT", "3306")),
    "user": os.getenv("MYSQL_USER", "root"),
    "password": os.getenv("MYSQL_PASSWORD", "12345678"),
    "charset": "utf8mb4",
    "cursorclass": pymysql.cursors.DictCursor,
}

QUESTION_IDS = [f"q{i}" for i in range(1, 11)]
METRIC_IDS = ["M1", "M2", "M3", "M4", "M5", "M6", "M7"]
QUESTIONS = set(QUESTION_IDS)
METRICS = set(METRIC_IDS)
AGENTS = ("agent1", "agent2", "agent3")
EXPECTED_ROW_COUNT = len(QUESTION_IDS) * len(METRIC_IDS)


def connect(database=True):
    cfg = dict(DB_CONFIG)
    if database:
        cfg["database"] = DB_NAME
    return pymysql.connect(**cfg)


def init_db():
    conn = connect(database=False)
    try:
        with conn.cursor() as c:
            c.execute(
                f"CREATE DATABASE IF NOT EXISTS `{DB_NAME}` "
                "CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci"
            )
        conn.commit()
    finally:
        conn.close()

    conn = connect(database=True)
    try:
        with conn.cursor() as c:
            c.execute(
                """
                CREATE TABLE IF NOT EXISTS scores (
                    id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
                    submission_id CHAR(36) NOT NULL,
                    scorer VARCHAR(100) NOT NULL,
                    role VARCHAR(100) NOT NULL DEFAULT '',
                    question VARCHAR(20) NOT NULL,
                    metric VARCHAR(10) NOT NULL,
                    agent1 TINYINT UNSIGNED NOT NULL,
                    agent2 TINYINT UNSIGNED NOT NULL,
                    agent3 TINYINT UNSIGNED NOT NULL,
                    comment TEXT NULL,
                    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
                    INDEX idx_question_metric (question, metric),
                    INDEX idx_metric (metric),
                    INDEX idx_submission (submission_id),
                    INDEX idx_created_at (created_at)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
                """
            )
            ensure_scores_schema(c)
        conn.commit()
    finally:
        conn.close()


def column_exists(cursor, table, column):
    cursor.execute(
        """
        SELECT COUNT(*) AS cnt
        FROM information_schema.COLUMNS
        WHERE TABLE_SCHEMA = %s AND TABLE_NAME = %s AND COLUMN_NAME = %s
        """,
        (DB_NAME, table, column),
    )
    return cursor.fetchone()["cnt"] > 0


def index_exists(cursor, table, index_name):
    cursor.execute(
        """
        SELECT COUNT(*) AS cnt
        FROM information_schema.STATISTICS
        WHERE TABLE_SCHEMA = %s AND TABLE_NAME = %s AND INDEX_NAME = %s
        """,
        (DB_NAME, table, index_name),
    )
    return cursor.fetchone()["cnt"] > 0


def ensure_scores_schema(cursor):
    if not column_exists(cursor, "scores", "submission_id"):
        cursor.execute(
            """
            ALTER TABLE scores
            ADD COLUMN submission_id CHAR(36) NOT NULL DEFAULT ''
            AFTER id
            """
        )

    if not column_exists(cursor, "scores", "question"):
        cursor.execute(
            """
            ALTER TABLE scores
            ADD COLUMN question VARCHAR(20) NOT NULL DEFAULT ''
            AFTER role
            """
        )

    if not column_exists(cursor, "scores", "comment"):
        cursor.execute(
            """
            ALTER TABLE scores
            ADD COLUMN comment TEXT NULL
            AFTER agent3
            """
        )

    if not column_exists(cursor, "scores", "updated_at"):
        cursor.execute(
            """
            ALTER TABLE scores
            ADD COLUMN updated_at TIMESTAMP NOT NULL
            DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
            AFTER created_at
            """
        )

    if not index_exists(cursor, "scores", "idx_question_metric"):
        cursor.execute("ALTER TABLE scores ADD INDEX idx_question_metric (question, metric)")

    if not index_exists(cursor, "scores", "idx_metric"):
        cursor.execute("ALTER TABLE scores ADD INDEX idx_metric (metric)")

    if not index_exists(cursor, "scores", "uq_scorer_role_question_metric"):
        try:
            cursor.execute(
                """
                ALTER TABLE scores
                ADD UNIQUE KEY uq_scorer_role_question_metric (scorer, role, question, metric)
                """
            )
        except pymysql.err.IntegrityError:
            # Existing duplicate historical rows should be resolved manually; do not delete user data here.
            pass


@app.route("/")
def index():
    return send_from_directory(DIR, "index.html")


@app.route("/images/<path:filename>")
def images(filename):
    return send_from_directory(os.path.join(DIR, "images"), filename)


@app.after_request
def cors(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS"
    return response


@app.route("/api/health")
def health():
    try:
        init_db()
        conn = connect()
        try:
            with conn.cursor() as c:
                c.execute("SELECT COUNT(*) AS total FROM scores")
                total = c.fetchone()["total"]
            return jsonify({"status": "ok", "database": DB_NAME, "scores": total})
        finally:
            conn.close()
    except Exception as e:
        return jsonify({"status": "error", "message": str(e)}), 500


def validate_score_row(row):
    question = str(row.get("question", "")).strip().lower()
    metric = str(row.get("metric", "")).strip().upper()
    if question not in QUESTIONS:
        raise ValueError(f"题目编号 {question or '<空>'} 无效")
    if metric not in METRICS:
        raise ValueError(f"{question} 的指标 {metric} 无效")

    values = {}
    for agent in AGENTS:
        try:
            value = int(row.get(agent))
        except (TypeError, ValueError):
            raise ValueError(f"{question} {metric} {agent} 必须是 0-10 的整数")
        if value < 0 or value > 10:
            raise ValueError(f"{question} {metric} {agent} 超出 0-10 范围")
        values[agent] = value

    return {
        "question": question,
        "metric": metric,
        "agent1": values["agent1"],
        "agent2": values["agent2"],
        "agent3": values["agent3"],
        "comment": str(row.get("comment", "")).strip(),
    }


def validate_complete_submission(rows):
    seen = set()
    for row in rows:
        key = (row["question"], row["metric"])
        if key in seen:
            raise ValueError(f"重复评分条目：{row['question'].upper()} {row['metric']}")
        seen.add(key)

    expected = {(question, metric) for question in QUESTION_IDS for metric in METRIC_IDS}
    missing = sorted(expected - seen, key=lambda item: (int(item[0][1:]), item[1]))
    extra = sorted(seen - expected, key=lambda item: (int(item[0][1:]), item[1]))
    if missing or extra or len(rows) != EXPECTED_ROW_COUNT:
        missing_text = "，".join(f"{q.upper()} {m}" for q, m in missing[:12])
        suffix = " ..." if len(missing) > 12 else ""
        raise ValueError(
            f"完整提交需要 {EXPECTED_ROW_COUNT} 条评分，当前 {len(rows)} 条。"
            f"缺少：{missing_text}{suffix}"
        )


@app.route("/api/scores", methods=["POST", "OPTIONS"])
def save_scores():
    if request.method == "OPTIONS":
        return ("", 204)

    data = request.get_json(silent=True) or {}
    scorer = str(data.get("scorer") or "").strip()[:100]
    role = str(data.get("role") or "").strip()[:100]
    comment = str(data.get("comment") or "").strip()
    scores = data.get("scores") or []

    if not scorer:
        return jsonify({"status": "error", "message": "请填写评估人姓名或编号"}), 400
    if not role:
        return jsonify({"status": "error", "message": "请选择评估人身份"}), 400

    if not isinstance(scores, list) or not scores:
        return jsonify({"status": "error", "message": "没有收到评分数据"}), 400

    try:
        rows = [validate_score_row(item) for item in scores]
        validate_complete_submission(rows)
    except ValueError as e:
        return jsonify({"status": "error", "message": str(e)}), 400

    init_db()
    submission_id = str(uuid.uuid4())
    conn = connect()
    try:
        with conn.cursor() as c:
            c.execute(
                """
                SELECT COUNT(*) AS existing
                FROM scores
                WHERE scorer = %s AND role = %s
                  AND question IN ('q1','q2','q3','q4','q5','q6','q7','q8','q9','q10')
                """,
                (scorer, role),
            )
            existing = c.fetchone()["existing"]
            if existing > 0:
                return (
                    jsonify(
                        {
                            "status": "error",
                            "message": "该评估人和身份已经提交过完整评分，不能重复提交。如需修改，请先由管理员清理该用户旧记录。",
                        }
                    ),
                    409,
                )
            for row in rows:
                row_comment = row["comment"] or comment
                c.execute(
                    """
                    INSERT INTO scores
                    (submission_id, scorer, role, question, metric, agent1, agent2, agent3, comment)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                    """,
                    (
                        submission_id,
                        scorer,
                        role,
                        row["question"],
                        row["metric"],
                        row["agent1"],
                        row["agent2"],
                        row["agent3"],
                        row_comment,
                    ),
                )
        conn.commit()
        return jsonify(
            {
                "status": "ok",
                "submission_id": submission_id,
                "count": len(rows),
                "mode": "created",
                "message": "已保存完整评分",
                "saved_at": datetime.now().isoformat(timespec="seconds"),
            }
        )
    except pymysql.err.IntegrityError:
        conn.rollback()
        return (
            jsonify(
                {
                    "status": "error",
                    "message": "该评估人和身份已经提交过评分，系统已阻止重复写入。",
                }
            ),
            409,
        )
    except Exception as e:
        conn.rollback()
        return jsonify({"status": "error", "message": str(e)}), 500
    finally:
        conn.close()


@app.route("/api/scores", methods=["GET"])
def get_scores():
    init_db()
    conn = connect()
    try:
        with conn.cursor() as c:
            c.execute(
                """
                SELECT question, metric,
                       AVG(agent1) AS avg_agent1,
                       AVG(agent2) AS avg_agent2,
                       AVG(agent3) AS avg_agent3,
                       COUNT(*) AS count
                FROM scores
                WHERE question IN ('q1','q2','q3','q4','q5','q6','q7','q8','q9','q10')
                GROUP BY question, metric
                ORDER BY CAST(SUBSTRING(question, 2) AS UNSIGNED), metric
                """
            )
            rows = c.fetchall()
        for row in rows:
            row["avg_agent1"] = round(float(row["avg_agent1"]), 2)
            row["avg_agent2"] = round(float(row["avg_agent2"]), 2)
            row["avg_agent3"] = round(float(row["avg_agent3"]), 2)
        return jsonify({"status": "ok", "data": rows})
    finally:
        conn.close()


@app.route("/api/results")
def results_page():
    init_db()
    conn = connect()
    try:
        with conn.cursor() as c:
            c.execute(
                """
                SELECT question, metric, AVG(agent1) AS a1, AVG(agent2) AS a2,
                       AVG(agent3) AS a3, COUNT(*) AS cnt
                FROM scores
                WHERE question IN ('q1','q2','q3','q4','q5','q6','q7','q8','q9','q10')
                GROUP BY question, metric
                ORDER BY CAST(SUBSTRING(question, 2) AS UNSIGNED), metric
                """
            )
            rows = c.fetchall()
    finally:
        conn.close()

    metric_names = {
        "M1": "回答相关性",
        "M2": "事实准确性",
        "M3": "信息完整性",
        "M4": "引用可核查性",
        "M5": "证据支撑程度",
        "M6": "幻觉控制",
        "M7": "科研实用性",
    }
    data = {}
    for row in rows:
        data.setdefault(row["question"], {})[row["metric"]] = (
            round(float(row["a1"]), 2),
            round(float(row["a2"]), 2),
            round(float(row["a3"]), 2),
            row["cnt"],
        )

    html = """
    <!DOCTYPE html><html lang="zh"><head><meta charset="UTF-8"><title>评估结果</title>
    <style>
    body{font:14px/1.6 'Segoe UI','Microsoft YaHei',sans-serif;max-width:960px;margin:36px auto;padding:0 20px;background:#f4f6f8;color:#263238}
    h1{color:#111827;border-bottom:2px solid #f59e0b;padding-bottom:10px}
    h2{margin-top:28px;color:#111827}
    table{width:100%;border-collapse:collapse;margin:12px 0 22px;background:#fff;box-shadow:0 1px 3px rgba(0,0,0,.08);border-radius:6px;overflow:hidden}
    td,th{border:1px solid #e5e7eb;padding:10px 12px;text-align:center}
    th{background:#111827;color:#fff;font-size:12px}
    td:first-child{text-align:left}
    .best{font-weight:700;color:#047857}
    .empty{padding:28px;background:#fff;border:1px solid #e5e7eb;border-radius:6px;color:#6b7280}
    </style></head><body><h1>Agent 人工评估结果</h1>
    """

    if not data:
        html += '<div class="empty">暂无评分数据。请先在首页提交评分。</div>'
    for question in QUESTION_IDS:
        if question not in data:
            continue
        html += f"<h2>{question.upper()}</h2><table><tr><th>指标</th><th>Agent 1</th><th>Agent 2</th><th>Agent 3</th><th>样本数</th></tr>"
        for metric in METRIC_IDS:
            if metric not in data[question]:
                continue
            a1, a2, a3, cnt = data[question][metric]
            best = max(a1, a2, a3)
            cells = []
            for value in (a1, a2, a3):
                cls = "best" if value == best else ""
                cells.append(f'<td class="{cls}">{value}</td>')
            html += f"<tr><td>{metric} {metric_names[metric]}</td>{''.join(cells)}<td>{cnt}</td></tr>"
        html += "</table>"
    html += "</body></html>"
    return html


if __name__ == "__main__":
    print("Server starting on http://localhost:5000")
    print("Health:      http://localhost:5000/api/health")
    print("Results:     http://localhost:5000/api/results")
    app.run(host="0.0.0.0", port=5000, debug=False)
