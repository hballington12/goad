from collections import OrderedDict
import os
from typing import List, Optional

from mkdocs_macros.plugin import MacrosPlugin

# Supported Languages and their metadata
LANGUAGES = OrderedDict(
    python={
        "extension": ".py",
        "display_name": "Python",
        "icon_name": "python",
        "code_name": "python",
    },
    rust={
        "extension": ".rs",
        "display_name": "Rust",
        "icon_name": "rust",
        "code_name": "rust",
    },
)


def code_tab(
    base_path: str,
    section: Optional[str],
    language_info: dict,
) -> str:
    """Generate a single tab for the code block corresponding to a specific language.
        It gets the code at base_path and possible section and pretty prints markdown for it

    Args:
        base_path (str): path where the code is located
        section (str, optional): section in the code that should be displayed
        language_info (dict): Language specific information (icon name, display name, ...)

    Returns:
        str: A markdown formatted string represented a single tab
    """
    language = language_info["code_name"]

    # Create path for Snippets extension
    snippets_file_name = f"{base_path}:{section}" if section else f"{base_path}"

    # See Content Tabs for details https://squidfunk.github.io/mkdocs-material/reference/content-tabs/
    return f"""=== \":fontawesome-brands-{language_info['icon_name']}: {language_info['display_name']}\"
    ```{language}
    --8<-- \"{snippets_file_name}\"
    ```
    """


def define_env(env: MacrosPlugin) -> None:
    @env.macro
    def code_block(
        path: str,
        section: str = None,  # pyright: ignore[reportArgumentType]
    ) -> str:
        """Dynamically generate a code block for the code located under {language}/path

        Args:
            path (str): base_path for each language
            section (str, optional): Optional segment within the code file. Defaults to None.
        Returns:
            str: Markdown tabbed code block
        """
        result = []

        for language, info in LANGUAGES.items():
            base_path = f"{language}/{path}{info['extension']}"
            full_path = "docs/source/src/" + base_path
            # Check if file exists for the language
            if os.path.exists(full_path):
                result.append(
                    code_tab(base_path, section, info)
                )

        return "\n".join(result)
