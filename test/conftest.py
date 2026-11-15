# SPDX-License-Identifier: BSD-3-Clause

import pooch
import pytest


@pytest.fixture(scope="session")
def registry():
    _registry = pooch.create(
        path=pooch.os_cache("osyris"),
        base_url="https://github.com/osyris-project/osyrisdata/releases/download/1.0.0/",  # noqa: E501
        registry={
            "cosmology.zip": "md5:a0ab245d9226e2d324291b27fa178021",
            "sedov.zip": "md5:363d4f1a61cf574638b9a0c39d1707e8",
            "sod.zip": "md5:6700035eb53f7848792919d3c6f9d772",
            "starformation.zip": "md5:3e92bc7389f3e5472157bf810230c0b3",
        },
    )
    return _registry


@pytest.fixture()
def cosmology_data_path(registry):
    return registry.fetch("cosmology.zip", processor=pooch.processors.Unzip())[0].split(
        "output_"
    )[0]


@pytest.fixture()
def sedov_data_path(registry):
    return registry.fetch("sedov.zip", processor=pooch.processors.Unzip())[0].split(
        "output_"
    )[0]


@pytest.fixture()
def sod_data_path(registry):
    return registry.fetch("sod.zip", processor=pooch.processors.Unzip())[0].split(
        "output_"
    )[0]


@pytest.fixture()
def starformation_data_path(registry):
    return registry.fetch("starformation.zip", processor=pooch.processors.Unzip())[
        0
    ].split("output_")[0]
