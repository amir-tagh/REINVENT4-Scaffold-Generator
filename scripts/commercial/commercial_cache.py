#!/usr/bin/env python3

"""Persistent SQLite cache for commercial molecular screening."""

import sqlite3
import hashlib
import json
import argparse

from pathlib import Path
from datetime import datetime, timezone

from molecular_identity import get_identity
from swisssimilarity_client import SwissSimilarityClient


DEFAULT_DB = (
    Path(__file__).resolve().parents[2]
    / "databases"
    / "commercial"
    / "screening_cache.sqlite"
)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


class CommercialCache:

    def __init__(self, database=DEFAULT_DB):

        self.database = Path(database)

        self.database.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.connection = sqlite3.connect(
            self.database,
            timeout=30,
        )

        self.connection.row_factory = sqlite3.Row

        self.connection.execute(
            "PRAGMA busy_timeout = 30000"
        )

        self.initialize()

    def initialize(self):

        with self.connection:

            self.connection.execute("""
                CREATE TABLE IF NOT EXISTS screening_cache (
                    cache_key TEXT PRIMARY KEY,
                    canonical_smiles TEXT NOT NULL,
                    inchikey TEXT NOT NULL,
                    library TEXT NOT NULL,
                    method TEXT NOT NULL,
                    status TEXT NOT NULL,
                    session_id TEXT,
                    raw_results TEXT,
                    hit_count INTEGER DEFAULT 0,
                    error_message TEXT,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
            """)

            self.connection.execute("""
                CREATE INDEX IF NOT EXISTS idx_inchikey
                ON screening_cache(inchikey)
            """)

            self.connection.execute("""
                CREATE INDEX IF NOT EXISTS idx_library
                ON screening_cache(library)
            """)

            self.connection.execute("""
                CREATE INDEX IF NOT EXISTS idx_status
                ON screening_cache(status)
            """)

    def close(self):
        self.connection.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    @staticmethod
    def make_key(smiles, library, method="ECFP"):

        identity = get_identity(smiles)

        payload = json.dumps(
            [
                "SwissSimilarity",
                identity.canonical_smiles,
                identity.inchikey,
                library,
                method,
            ],
            separators=(",", ":"),
        )

        key = hashlib.sha256(
            payload.encode("utf-8")
        ).hexdigest()

        return key, identity

    def lookup(self, smiles, library, method="ECFP"):

        key, _ = self.make_key(
            smiles,
            library,
            method,
        )

        row = self.connection.execute("""
            SELECT *
            FROM screening_cache
            WHERE cache_key = ?
        """, (key,)).fetchone()

        return dict(row) if row else None

    def save_session(
        self,
        smiles,
        library,
        method,
        session_id,
    ):

        key, identity = self.make_key(
            smiles,
            library,
            method,
        )

        now = utc_now()

        with self.connection:

            self.connection.execute("""
                INSERT INTO screening_cache (
                    cache_key,
                    canonical_smiles,
                    inchikey,
                    library,
                    method,
                    status,
                    session_id,
                    created_at,
                    updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)

                ON CONFLICT(cache_key) DO UPDATE SET
                    status = excluded.status,
                    session_id = excluded.session_id,
                    raw_results = NULL,
                    hit_count = 0,
                    error_message = NULL,
                    updated_at = excluded.updated_at
            """, (
                key,
                identity.canonical_smiles,
                identity.inchikey,
                library,
                method,
                "SUBMITTED",
                str(session_id),
                now,
                now,
            ))

        return key

    def save_results(
        self,
        smiles,
        library,
        method,
        session_id,
        raw_results,
    ):

        key, identity = self.make_key(
            smiles,
            library,
            method,
        )

        hits = SwissSimilarityClient.parse_results(
            raw_results,
            library,
        )

        now = utc_now()

        with self.connection:

            self.connection.execute("""
                INSERT INTO screening_cache (
                    cache_key,
                    canonical_smiles,
                    inchikey,
                    library,
                    method,
                    status,
                    session_id,
                    raw_results,
                    hit_count,
                    created_at,
                    updated_at
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)

                ON CONFLICT(cache_key) DO UPDATE SET
                    status = excluded.status,
                    session_id = excluded.session_id,
                    raw_results = excluded.raw_results,
                    hit_count = excluded.hit_count,
                    error_message = NULL,
                    updated_at = excluded.updated_at
            """, (
                key,
                identity.canonical_smiles,
                identity.inchikey,
                library,
                method,
                "COMPLETED",
                str(session_id),
                raw_results,
                len(hits),
                now,
                now,
            ))

        return key

    def save_error(
        self,
        smiles,
        library,
        method,
        error_message,
    ):

        key, _ = self.make_key(
            smiles,
            library,
            method,
        )

        with self.connection:

            cursor = self.connection.execute("""
                UPDATE screening_cache
                SET status = ?,
                    error_message = ?,
                    updated_at = ?
                WHERE cache_key = ?
                  AND status != 'COMPLETED'
            """, (
                "ERROR",
                str(error_message),
                utc_now(),
                key,
            ))

        return cursor.rowcount > 0

    def get_completed(
        self,
        smiles,
        library,
        method="ECFP",
    ):

        record = self.lookup(
            smiles,
            library,
            method,
        )

        if record is None:
            return None

        if record["status"] != "COMPLETED":
            return None

        return record

    def statistics(self):

        total = self.connection.execute("""
            SELECT COUNT(*) FROM screening_cache
        """).fetchone()[0]

        molecules = self.connection.execute("""
            SELECT COUNT(DISTINCT canonical_smiles)
            FROM screening_cache
        """).fetchone()[0]

        hits = self.connection.execute("""
            SELECT COALESCE(SUM(hit_count), 0)
            FROM screening_cache
            WHERE status = 'COMPLETED'
        """).fetchone()[0]

        statuses = self.connection.execute("""
            SELECT status, COUNT(*) AS count
            FROM screening_cache
            GROUP BY status
        """).fetchall()

        libraries = self.connection.execute("""
            SELECT library, COUNT(*) AS count
            FROM screening_cache
            GROUP BY library
        """).fetchall()

        return {
            "Database": str(self.database),
            "Total_Records": total,
            "Unique_Query_Molecules": molecules,
            "Total_Retrieved_Hits": hits,
            "Status_Counts": {
                row["status"]: row["count"]
                for row in statuses
            },
            "Library_Counts": {
                row["library"]: row["count"]
                for row in libraries
            },
        }


def main():

    parser = argparse.ArgumentParser(
        description="Commercial screening SQLite cache"
    )

    parser.add_argument(
        "--database",
        default=str(DEFAULT_DB),
    )

    parser.add_argument(
        "--stats",
        action="store_true",
    )

    args = parser.parse_args()

    with CommercialCache(args.database) as cache:

        if args.stats:

            print(json.dumps(
                cache.statistics(),
                indent=4,
            ))

        else:

            print(
                "Cache initialized:",
                cache.database,
            )


if __name__ == "__main__":
    main()
