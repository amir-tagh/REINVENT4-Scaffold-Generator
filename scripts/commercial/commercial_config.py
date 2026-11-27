#!/usr/bin/env python3

"""Configuration management for commercial screening."""

import argparse
import json
from pathlib import Path
from collections import Counter


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_CONFIG = (
    PROJECT_ROOT / "config" / "commercial_sources.json"
)


class CommercialConfig:

    def __init__(self, filename=DEFAULT_CONFIG):

        self.filename = Path(filename)

        with self.filename.open() as handle:
            self.data = json.load(handle)

        self.validate()

    def validate(self):

        required = {
            "schema_version",
            "provider",
            "base_url",
            "screening",
            "api",
            "cache",
            "libraries",
            "profiles",
            "default_profile",
        }

        missing = required - set(self.data)

        if missing:
            raise ValueError(
                f"Missing configuration sections: {missing}"
            )

        libraries = self.data["libraries"]
        profiles = self.data["profiles"]

        for name, members in profiles.items():

            if not members:
                raise ValueError(
                    f"Empty profile: {name}"
                )

            if len(members) != len(set(members)):
                raise ValueError(
                    f"Duplicate libraries in profile: {name}"
                )

            for library in members:

                if library not in libraries:
                    raise ValueError(
                        f"Unknown library {library} "
                        f"in profile {name}"
                    )

        if self.data["default_profile"] not in profiles:
            raise ValueError("Invalid default profile")

        threshold = self.data["screening"][
            "minimum_analog_similarity"
        ]

        if not 0 <= threshold <= 1:
            raise ValueError("Invalid similarity threshold")

        if self.data["screening"]["mode"] not in {
            "exact", "similarity", "both"
        }:
            raise ValueError("Invalid screening mode")

        if self.data["api"]["parallel_jobs"] != 1:
            raise ValueError(
                "Only sequential screening is supported"
            )

        for option in (
            "request_timeout",
            "poll_interval",
            "maximum_wait",
            "request_interval",
        ):
            if self.data["api"][option] <= 0:
                raise ValueError(
                    f"Invalid API setting: {option}"
                )

        if self.data["api"]["maximum_new_jobs_per_run"] < 0:
            raise ValueError("Invalid new-job limit")

        if self.data["screening"][
            "maximum_analogs_per_library"
        ] < 0:
            raise ValueError("Invalid analog limit")

        if not self.data["base_url"].startswith("https://"):
            raise ValueError("HTTPS endpoint required")

    @property
    def libraries(self):
        return self.data["libraries"]

    @property
    def screening(self):
        return self.data["screening"]

    @property
    def api(self):
        return self.data["api"]

    @property
    def cache(self):
        return self.data["cache"]

    @property
    def base_url(self):
        return self.data["base_url"]

    def get_libraries(
        self,
        profile=None,
        enabled_only=True,
    ):

        profile = profile or self.data["default_profile"]

        if profile == "all":
            selected = list(self.libraries)

        else:

            if profile not in self.data["profiles"]:
                raise ValueError(
                    f"Unknown screening profile: {profile}"
                )

            selected = self.data["profiles"][profile]

        if enabled_only:

            selected = [
                library
                for library in selected
                if self.libraries[library]["enabled"]
            ]

        return selected

    def database_path(self):

        path = Path(self.cache["database"])

        if not path.is_absolute():
            path = PROJECT_ROOT / path

        return path

    def summary(self):

        categories = Counter(
            item["category"]
            for item in self.libraries.values()
        )

        return {
            "Provider": self.data["provider"],
            "Total_Libraries": len(self.libraries),
            "Enabled_Libraries": [
                key
                for key, value in self.libraries.items()
                if value["enabled"]
            ],
            "Validated_Libraries": [
                key
                for key, value in self.libraries.items()
                if value["validated"]
            ],
            "Categories": dict(categories),
            "Default_Profile": self.data["default_profile"],
            "Screening_Method": self.screening["method"],
            "Screening_Mode": self.screening["mode"],
            "Minimum_Analog_Similarity":
                self.screening["minimum_analog_similarity"],
            "Remote_Submission_Allowed":
                self.api["allow_remote_submission"],
            "Cache_Database": str(self.database_path()),
        }


def main():

    parser = argparse.ArgumentParser(
        description="Commercial screening configuration"
    )

    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG),
    )

    parser.add_argument(
        "--profile",
        default=None,
    )

    parser.add_argument(
        "--include-disabled",
        action="store_true",
    )

    args = parser.parse_args()

    config = CommercialConfig(args.config)

    print(json.dumps(
        config.summary(),
        indent=4,
    ))

    print("\nSelected libraries:")

    for library in config.get_libraries(
        profile=args.profile,
        enabled_only=not args.include_disabled,
    ):
        print(" ", library)


if __name__ == "__main__":
    main()
