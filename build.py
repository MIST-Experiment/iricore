from typing import Any, Dict

from setuptools_cpp import CMakeExtension, ExtensionBuilder

ext_modules = [
    CMakeExtension(f"libiri", sourcedir="src/iricore"),
]


def build(setup_kwargs: Dict[str, Any]) -> None:
    from src.iricore.data_update import update_iri
    update_iri()
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmdclass": dict(build_ext=ExtensionBuilder),
            "zip_safe": False,
        }
    )
