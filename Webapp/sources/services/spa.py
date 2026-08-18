import json
import os

from flask import current_app, url_for
from markupsafe import Markup


def vite_entry(entry):
    if current_app.config.get("VITE_DEV_SERVER", False):
        return _vite_dev_entry(entry)

    return _vite_production_entry(entry)


def _vite_dev_entry(entry):
    """
    Injects details of local vite server for quick refresh
    Args:
        entry:

    Returns:

    """
    vite_url = current_app.config.get(
        "VITE_DEV_SERVER_URL",
        "http://localhost:8000",
    ).rstrip("/")

    return Markup(  # noqa
        f"""
<script type="module">
    import RefreshRuntime from "{vite_url}/@react-refresh";

    RefreshRuntime.injectIntoGlobalHook(window);
    window.$RefreshReg$ = () => {{}};
    window.$RefreshSig$ = () => (type) => type;
    window.__vite_plugin_react_preamble_installed__ = true;
</script>

<script
    type="module"
    src="{vite_url}/@vite/client"
></script>

<script
    type="module"
    src="{vite_url}/{entry}"
></script>
"""
    )


def _vite_production_entry(entry):
    manifest_path = os.path.join(
        current_app.root_path,
        "static",
        "spa",
        ".vite",
        "manifest.json",
    )

    with open(manifest_path, encoding="utf-8") as manifest_file:
        manifest = json.load(manifest_file)

    asset = manifest[entry]
    imported_chunks = _get_imported_chunks(manifest, entry)

    tags = []
    seen_css = set()

    # Entry CSS
    for css_file in asset.get("css", []):
        if css_file not in seen_css:
            seen_css.add(css_file)

            tags.append(
                '<link rel="stylesheet" href="{}">'.format(
                    url_for(
                        "static",
                        filename=f"spa/{css_file}",
                    )
                )
            )

    # CSS belonging to imported chunks
    for chunk in imported_chunks:
        for css_file in chunk.get("css", []):
            if css_file not in seen_css:
                seen_css.add(css_file)

                tags.append(
                    '<link rel="stylesheet" href="{}">'.format(
                        url_for(
                            "static",
                            filename=f"spa/{css_file}",
                        )
                    )
                )

    # Main entry
    tags.append(
        '<script type="module" src="{}"></script>'.format(
            url_for(
                "static",
                filename=f"spa/{asset['file']}",
            )
        )
    )

    return Markup("\n".join(tags))


def _get_imported_chunks(manifest, entry):
    seen = set()
    chunks = []

    def collect(chunk):
        for import_name in chunk.get("imports", []):
            if import_name in seen:
                continue

            seen.add(import_name)

            imported_chunk = manifest[import_name]

            collect(imported_chunk)
            chunks.append(imported_chunk)

    collect(manifest[entry])

    return chunks
