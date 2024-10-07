"""
Parse XML tags for easy styling of CLI text output.
---------------------------------------------------
Author: Moenen Erbuer - moenen.erbuer@ibm.com
v0.0.0-beta6 / Last update: Sep 28, 2023

Description:
    This module parses XML style tags into ANSI escape codes,
    allowing for easy styling of CLI text output (color, bold, etc.)

Available functions:
    style()                 Returns styled text
    print_s()               Print styled text
    strip_tags()            Remove style tags
    tags_to_markdown()      Convert tags to iPython Markdown (for use in Jupyter)
    wrap_text()             Limit paragraph width

Basic usage:
    print_s('<error>Something's wrong</error>', tab=1, edge=True, pad=2)
    foo = style('<error>Something's wrong</error>', tab=1, edge=True, pad=2)
    print(foo)

Note:
    The <cmd> tag is used for displaying commands.
    Command parameters are displayed in <angle_brackets>,
    so non-recognized tags that appear inside <cmd> tags
    are displayed soft:
    <cmd>open file <file_name></cmd>

To do:
    Update wrap_text to account for ANSI escape codes - see function in JSON editor.

"""

import sys
import re
import shutil
from collections import deque

# Style tags.
# - - -
# Any tag that is not in this dictionary will be ignored.
# When a non-recognized tag is found inside the <cmd> tag,
# it is displayed as <soft>, because that's how we notate
# parameters in our help functionality.
tags = {
    # Custom tags
    "h1": "\x1b[0m",  # Primary headers: default with yellow line below
    "h2": "\x1b[33m",  # Secondary headers: yellow
    "cmd": "\x1b[36m",  # Commands: cyan
    "error": "\x1b[31m",  # Errors: red
    "error_reverse": "\x1b[37;41m",  # Error reverse: white on red
    "warning": "\x1b[33m",  # Warnings: yellow
    "success": "\x1b[32m",  # Success: green
    "link": "\x1b[4;94m",  # Links: undertline intense blue
    "edit": "\x1b[7m",  # Edit node: reverse
    "editable": "\x1b[0;47;30m",  # Editable text: black + on_white
    # Styles
    "reset": "\x1b[0m",
    "bold": "\x1b[1m",
    "dim": "\x1b[2m",  # Current color but dimmed
    "soft": "\x1b[90m",  # Gray
    "italic": "\x1b[3m",
    "underline": "\x1b[4m",
    "blink": "\x1b[5m",
    "reverse": "\x1b[7m",
    "hidden": "\x1b[8m",
    "strikethrough": "\x1b[9m",
    # Foreground colors - Regular
    "black": "\x1b[30m",
    "red": "\x1b[31m",
    "green": "\x1b[32m",
    "yellow": "\x1b[33m",
    "blue": "\x1b[34m",
    "magenta": "\x1b[35m",
    "cyan": "\x1b[36m",
    "white": "\x1b[37m",
    # Foreground colors - Bright/intense (non-standard)
    "bright_black": "\x1b[90m",
    "bright_red": "\x1b[91m",
    "bright_green": "\x1b[92m",
    "bright_yellow": "\x1b[93m",
    "bright_blue": "\x1b[94m",
    "bright_magenta": "\x1b[95m",
    "bright_cyan": "\x1b[96m",
    "bright_white": "\x1b[97m",
    # Background colors - Regular
    "on_black": "\x1b[0;40m",
    "on_red": "\x1b[0;41m",
    "on_green": "\x1b[0;42m",
    "on_yellow": "\x1b[0;43m",
    "on_blue": "\x1b[0;44m",
    "on_magenta": "\x1b[0;45m",
    "on_cyan": "\x1b[0;46m",
    "on_white": "\x1b[0;47m",
    # Background colors - Bright/intense (non-standard)
    "on_bright_black": "\x1b[100m",
    "on_bright_red": "\x1b[101m",
    "on_bright_green": "\x1b[102m",
    "on_bright_yellow": "\x1b[103m",
    "on_bright_blue": "\x1b[104m",
    "on_bright_magenta": "\x1b[105m",
    "on_bright_cyan": "\x1b[106m",
    "on_bright_white": "\x1b[107m",
}


def style(
    text: str,
    pad: int = 0,
    pad_top: int = 0,
    pad_btm: int = 0,
    tabs: int = 0,
    edge: bool | str = False,
    width: int = None,
    nowrap: bool = False,
    xml_tags: bool = True,
    trim: bool = False,
):
    """
    Parse xml tags and return ANSI infused output.

    Parameters:
        • text (str): Input text.
        • pad (int): Number of line breaks before and after.
        • pad_top (int): Number of line breaks before.
        • pad_btm (int): Number of line breaks after.
        • tabs (int): Number of tabs to indent.
        • edge (bool|str): Line left edge with | Can be color or True for grey.
        • width (int): Width of paragraph in characters.
        • nowrap (bool): Do not wrap lines.
        • xml_tags (bool): Whether to convert XML tags.
        • trim (bool): Remove line breaks at front/end.

    Returns:
        Text encoded with ANSI escapse codes ready to be printed to the CLI.
    """
    if width is None:
        try:
            columns = shutil.get_terminal_size().columns
            width = min(columns - 5, 150)
            if tabs:
                width = width - (tabs * 4)
            elif edge:
                width = width - 5
        except Exception:  # pylint: disable=broad-except
            width = 60

    if text is None:
        return ""

    # Enforce string type.
    text = str(text)

    # Strip any line breaks at start/end.
    if trim:
        text = _trim(text)

    # Add dotted lines under h1 tags.
    text = _add_header_lines(text)

    # Replace line breaks so all text is parsed on one line.
    text = text.replace("\n", "---LINEBREAK1---")

    # Make links clickable.
    text = _parse_links(text)

    if xml_tags == True:
        # Style xml command parameters.
        text = _style_cmd_params(text)

        # Style xml style tags.
        text = _parse_tags(text)

    # Restore line breaks.
    text = text.replace("---LINEBREAK1---", "\n")

    # Add tabs to the start of each line.
    if tabs or edge:
        tabs = tabs if tabs else 1
        text = _add_tabs(text, tabs)

    # Wrap lines based on paragraph width in characters.
    if not nowrap:
        text = wrap_text(text, width=width)

    # Add edge.
    if edge:
        text = _edge(text, edge)

    # Add top and bottom padding
    if pad and pad != 0 and isinstance(pad, int):
        padding = "\n" * pad
        text = padding + text + padding
    else:
        if pad_top:
            padding = "\n" * pad_top
            text = padding + text
        if pad_btm:
            padding = "\n" * pad_btm
            text = text + padding

    return text


def print_s(text: str, ephemeral=False, pre_styled_bulk=False, **kwargs):
    """
    Print styled text. This is a wrapper around style().

    Parameters:
        • See style()
        • ephemeral (bool): Print next line on the same line. (Currently not used, may act wonky)
    """

    text = str(text)
    end = "\r" if ephemeral else "\n"
    if pre_styled_bulk is True:
        sys.stdout.write(text)
        sys.stdout.flush()
        print("")
    else:
        print(style(text, **kwargs), end=end, flush=True)


def strip_tags(text: str):
    """Recursively remove all XML tags."""

    def _strip(match: object, pattern):
        inner_text = match.group(2)

        # Strip any nested tags.
        if re.findall(pattern, inner_text):
            inner_text = re.sub(pattern, lambda match: _strip(match, pattern), inner_text)

        return inner_text

    #
    #

    # Enforce string type.
    text = str(text)

    # Replace line breaks so all text is parsed on one line.
    text = text.replace("\n", "---LINEBREAK2---")

    # Expand error and success tags.
    text = _expand_error_success_tags(text, False)

    # Replace recognized tags with appropriate styling.
    pattern = rf"<({'|'.join(list(tags))})>(.*?)</\1>"
    text = re.sub(pattern, lambda match: _strip(match, pattern), text)

    # Restore line breaks.
    text = text.replace("---LINEBREAK2---", "\n")

    return text


# UNTESTED AND UNUSED
# def strip_ansi(text: str):
#     """
#     Remove ANSI escape codes.

#     Code from https://stackoverflow.com/a/2187024
#     """

#     from pyparsing import Literal, Word, Combine, Optional, delimitedList, oneOf, alphas, Suppress, nums

#     ESC = Literal('\x1b')
#     integer = Word(nums)
#     escapeSeq = Combine(ESC + '[' + Optional(delimitedList(integer, ';')) + oneOf(list(alphas)))

#     def nonAnsiString(s): return Suppress(escapeSeq).transformString(s)

#     unColorString = nonAnsiString('\x1b[1m0.0\x1b[0m')
#     print(unColorString, len(unColorString))


def tags_to_markdown(text: str):
    if text is None:
        return ""

    # Enforce string type.
    text = str(text)

    # Replace leading spaces with non-breaking spaces.
    # text = re.sub(r'^( *)', '', text, flags=re.MULTILINE) # Forgot why we needed this, but it's fucking with ASCII type in some splash pages (%openad openad)

    # Replace line breaks so all text is parsed on one line.
    # Because html breaks (<br>) don't play well with headings,
    # and end of line characters don't play well with `code`
    # blocks, we have to do some trickery here.
    text = re.sub(r"(</h[123]>)(\n+)", lambda match: match.group(1) + len(match.group(2)) * "---LINEBREAKSOFT---", text)
    text = re.sub(r"(\n+)(<h[123]>)", lambda match: len(match.group(1)) * "---LINEBREAKSOFT---" + match.group(2), text)
    text = _replace_linebreaks_inside_cmdblocks(text, "---LINEBREAK3---")
    # Every linebreak becomes a hard <br> except when the line has no text.
    text = "".join([line + "---LINEBREAK3---" if line.strip() else "---LINEBREAKSOFT---" for line in text.splitlines()])

    # Expand error and success tags.
    text = _expand_error_success_tags(text, True)

    # Replace <soft> and <underline> tags only
    # when they don't appear inside <cmd> tags.
    text = _replace_soft_and_underline(text)

    # Replace all other tags.
    text = re.sub(r"<h1>(.*?)<\/h1>", r"## \1", text)
    text = re.sub(r"<h2>(.*?)<\/h2>", r"### \1", text)
    text = re.sub(r"<link>(.*?)<\/link>", r'<a target="_blank" href="\1">\1</a>', text)
    text = re.sub(r"<bold>(.*?)<\/bold>", r"**\1**", text)
    text = re.sub(r"<cmd>(.*?)<\/cmd>", r"`\1`", text)
    text = re.sub(r"<red>(.*?)<\/red>", r'<span style="color: #d00">\1</span>', text)
    text = re.sub(r"<green>(.*?)<\/green>", r'<span style="color: #090">\1</span>', text)
    text = re.sub(r"<yellow>(.*?)<\/yellow>", r'<span style="color: #dc0">\1</span>', text)
    text = re.sub(r"<blue>(.*?)<\/blue>", r'<span style="color: #00d">\1</span>', text)
    text = re.sub(r"<magenta>(.*?)<\/magenta>", r'<span style="color: #d07">\1</span>', text)
    text = re.sub(r"<cyan>(.*?)<\/cyan>", r'<span style="color: #0cc">\1</span>', text)
    text = re.sub(r"<on_red>(.*?)<\/on_red>", r'<span style="background: #d00; color: #fff">\1</span>', text)
    text = re.sub(r"<on_green>(.*?)<\/on_green>", r'<span style="background: #090; color: #fff">\1</span>', text)
    text = re.sub(r"<on_yellow>(.*?)<\/on_yellow>", r'<span style="background: #dc0; color: #fff">\1</span>', text)
    text = re.sub(r"<on_blue>(.*?)<\/on_blue>", r'<span style="background: #00d; color: #fff">\1</span>', text)
    text = re.sub(r"<on_magenta>(.*?)<\/on_magenta>", r'<span style="background: #d07; color: #fff">\1</span>', text)
    text = re.sub(r"<on_cyan>(.*?)<\/on_cyan>", r'<span style="background: #0cc; color: #fff">\1</span>', text)

    # Remove empty code tags.
    # These can be generated when having style tags
    # (which will be removed) inside of a command tag.
    # For example "<cmd>foo <soft>bar</soft></cmd>" will
    # become "`foo `bar``" which can cause problems
    # when it's followed by a linebreak.
    text = re.sub(r"``", "", text)

    # Escape quotes.
    text = text.replace("'", "'")

    # Remove all other tags
    text = strip_tags(text)

    # Restore line breaks.
    text = text.replace("---LINEBREAKSOFT---", "\n")
    text = text.replace("---LINEBREAK3---", " <br> \n")

    return text


# In Jupyter, when you render codeblocks using `` quotes,
# the <br> linebreaks are rendered as raw text. To fix this,
# we replace `foo<br>bar` with `foo`<br>`bar`.
# This function is run before we replace \n with <br>, so
# we replace \n instead of <br>
def _replace_linebreaks_inside_cmdblocks(text: str, break_str: str):
    # First strip outer line breaks, to avoid parsing these:
    # <cmd>
    #   foo
    # </cmd>
    pattern1 = r"<cmd>(\n)?([^<]*?)(\n)?</cmd>"
    text = re.sub(pattern1, r"<cmd>\2</cmd>", text)

    # Next replace line breaks inside <cmd> tags.
    pattern2 = r"<cmd>([^<]*?)\n([^<]*?)</cmd>"
    while re.search(pattern2, text):
        text = re.sub(pattern2, rf"<cmd>\1</cmd>{break_str}<cmd>\2</cmd>", text)
    return text


# Replace <soft> and <underline> tags only
# when they don't appear inside <cmd> tags.
def _replace_soft_and_underline(text: str):
    # Remove the tags within <cmd> tags.
    text = re.sub(
        r"<cmd>(.*?)<\/cmd>", lambda x: x.group(0).replace("<underline>", "").replace("</underline>", ""), text
    )
    text = re.sub(r"<cmd>(.*?)<\/cmd>", lambda x: x.group(0).replace("<soft>", "`").replace("</soft>", "`"), text)

    # Replace the tags outside of <cmd> tags.
    text = re.sub(r"<underline>(.*?)<\/underline>", r'<span style="text-decoration: underline">\1</span>', text)
    text = re.sub(r"<soft>(.*?)<\/soft>", r'<span style="color: #ccc">\1</span>', text)
    return text


def wrap_text(text: str, width=80):
    """
    Wrap lines to fit within the paragraph width.

    Parameters:
        • text (str): Input text.
        • width (int): Paragraph width in characters.
    """

    text = str(text)
    lines = text.splitlines()
    wrapped_lines = []

    for line in lines:
        stripped_line = line.lstrip()
        indent = line[: len(line) - len(stripped_line)]
        wrapped = a_textwrap(stripped_line, width=width, subsequent_indent=indent)
        wrapped_lines.append(indent + wrapped)

    return "\n".join(wrapped_lines)


def a_len(text):
    """An alternative to len() which measures a string's 'actual' length without ANSI codes."""
    pattern_ansi = r"\x1b\[[0-9;]*[m]"
    return len(re.sub(pattern_ansi, "", text))


def a_textwrap(text, width=100, drop_whitespace=True, subsequent_indent=""):
    """An alternative to textwrap.fill() that accounts for ANSI codes."""

    # Split by ANSI characters and spaces.
    # By using capturing groups, the separators are included in the list.
    pattern_ansi_or_space_captured = r"(\x1b\[[0-9;]*[m])|([\s]+)"
    split_text = re.split(pattern_ansi_or_space_captured, text)
    split_text = list(filter(None, split_text))  # Remove the empty strings.

    output = []
    line = ""  # With ansi codes, used for output.

    def _linebreak():
        nonlocal output, line
        output.append(line)
        line = subsequent_indent

    for i, string in enumerate(split_text):
        pattern_ansi = r"\x1b\[[0-9;]*[m]"
        is_ansi = re.match(pattern_ansi, string)
        str_len_clean = len(re.sub(pattern_ansi, "", string))
        reset_ansi_code = "\x1b[0m"

        if string == "\n":
            # Add linebreak.
            _linebreak()
        elif is_ansi:
            # Add ANSI code to line, add linebreak before if the line is full.
            if a_len(line) == width and string != reset_ansi_code:
                _linebreak()
            line += string
        elif a_len(line) + str_len_clean <= width:
            # Add string to line when there's room,
            # but only if the string is not a space,
            # or a space under certain circumstances.
            line += string
        else:
            # Start new line.
            _linebreak()
            if string and (not string.isspace() or not drop_whitespace):
                line += string

        # Wrap up last line.
        if i == len(split_text) - 1:
            output.append(line)

    return "\n".join(output)


#
#
#
#


def _trim(text: str):
    """Remove line breaks at front/end caused by putting quotes on separate line."""
    lines = deque(text.split("\n"))

    if lines and lines[0].strip() == "":
        lines.popleft()
    if lines and lines[len(lines) - 1].strip() == "":
        lines.pop()

    return "\n".join(lines)


def _add_header_lines(text):
    """Add dotted lines under h1 tags."""
    pattern = rf"(.*?)<h1>(.*?)</h1>"

    def _modify(match: object):
        original_text = match.group(0)
        space_before = match.group(1)
        inner_text = match.group(2)
        lines = space_before + "<yellow>" + ("-" * len(strip_tags(inner_text))) + "</yellow>"
        return original_text + "\n" + lines

    return re.sub(pattern, _modify, text)


def _parse_links(text):
    """Parse html links and return clickable output."""
    lines = text.splitlines()
    result = []
    for line in lines:
        # Detect <a href="x.yz">click here</a>
        pattern = r'<a href="(.*)">(.*?)<\/a>'
        link_tuples = re.findall(pattern, line)  # [(url, text), ...]
        # ['<a href="xyz">click here</a>', ...]
        match_strings = re.search(pattern, line)
        if link_tuples:
            for item in enumerate(link_tuples):
                i = item[0]
                link_tuple = item[1]

                # Create clickable URL.
                url, text = link_tuple
                link_interactive = f"\033]8;;{url}\033\\{text}\033]8;;\033\\"

                # Replace link string with clickable URL.
                link_string = match_strings[i]
                result.append(line.replace(link_string, link_interactive))
        else:
            result.append(line)

    return "\n".join(result)


def _style_cmd_params(text: str):
    """
    Mark non-style related tags as soft.

    Some context:
    We use XML tags to style text: <red>red text</red>
    These tags get replaced with ANSI escape codes in parse_tags().
    However, command parameters are noted in the same angle bracket
    format, eg. `run <run_name>`, and displayed as soft text, making
    it easier to read the commands.

    Because they always appear inside the <cmd> tag, we first need to close
    the cmd tag, then open the soft tag, then open the cmd tag again.

    So in practice:
    <cmd>run <run_name></cmd> --> <cmd>run </cmd><soft><run_name></soft><cmd>
    """

    text = str(text)
    lines = text.splitlines()
    result = []

    for line in lines:
        # Find all tags.
        line_copy = line
        searching = True
        params = []
        while searching:
            pattern = rf'<cmd>.*?(<(?!{"|".join(list(tags))})[^</>]+>).*</cmd>'
            param = re.findall(pattern, line_copy)
            if param:
                param = param[0]
                params.append(param)
                line_copy = line_copy.replace(param, "")
                line = line.replace(param, f"<soft>{param}</soft>")
            else:
                searching = False

        result.append(line)

    return "\n".join(result)


def _parse_tags(text: str):
    """Parse xml tags and return styled output."""

    # Enforce string type.
    text = str(text)

    # Remove <pre> tags which are used to wrap ASCII headers (see splash.py)
    text = re.sub(r"<pre[^<]*?>(.*?)</pre>", r"\1", text)  # flags=re.MULTILINE

    # Replace open-close tags <red>lorem ipsum</red>.
    pattern = rf"(.*?)<({'|'.join(list(tags))})>(.*?)</\2>"
    text = re.sub(pattern, lambda match: __replace(match, pattern), text)

    # Replace singular tags <red/> --> everything red from here on.
    pattern_singular = r"<({})/>".format("|".join(list(tags)))
    text = re.sub(pattern_singular, __replace_singular, text)
    return text


def __replace(match: object, pattern, parent_tag="reset"):
    """Replace regex matches with appropriate styling."""
    text_before = match.group(1)
    tag = match.group(2)
    inner_text = match.group(3)
    ansi_code_open = tags[tag]
    ansi_code_close = tags[parent_tag]
    # We need ansi_code_reset before ansi_code_close
    # for underline/bold to work.
    ansi_code_reset = tags["reset"]

    # Replace any nested tags.
    if re.findall(pattern, inner_text):
        inner_text = re.sub(pattern, lambda match: __replace(match, pattern, tag), inner_text)

    return f"{text_before}{ansi_code_open}{inner_text}{ansi_code_reset}{ansi_code_close}"


def __replace_singular(match: object):
    tag = match.group(1)
    return tags[tag]


def _add_tabs(text: str, tabs: int):
    """
    Add left padding to a text block.

    Parameters:
        • text (str): Input text.
        • tabs (int): Number of tabs to offset, each tab being 4 spaces.
    """

    text = str(text)
    spacing = "    " * tabs
    return spacing + ("\n" + spacing).join(text.split("\n"))


def _edge(text, edge):
    """Add edge to left side of a text block."""
    edge_color = "soft" if edge == True else edge
    edge_str = "|"
    lines = text.splitlines()
    placc = "\x1b[0m"  # Previous line ANSI color character
    for i, line in enumerate(lines):
        edge_str = tags[edge_color] + edge_str + placc
        placc = re.findall(r"\x1b\[\d+m", line)
        placc = placc.pop() if len(placc) else ""
        lines[i] = edge_str + line

    return "\n".join(lines)


def _expand_error_success_tags(text, html=False):
    """Turn error/success tags into text, for non-styled output."""

    # Output for markdown: use HTML tags instead of ANSI codes.
    if html:
        text = re.sub(r"<error>(.*?)<\/error>", r'<span style="color: #d00">\1</span>', text, flags=re.MULTILINE)
        text = re.sub(r"<success>(.*?)<\/success>", r'<span style="color: #090">\1</span>', text, flags=re.MULTILINE)
        text = re.sub(r"<warning>(.*?)<\/warning>", r'<span style="color: #ffa500">\1</span>', text, flags=re.MULTILINE)

    # Output for non-styled output: display status as text.
    else:
        text = re.sub(r"^<error>(.*?)<\/error>", r"Error: \1", text, flags=re.MULTILINE)
        text = re.sub(r"^<success>(.*?)<\/success>", r"Success: \1", text, flags=re.MULTILINE)
        text = re.sub(r"^<warning>(.*?)<\/warning>", r"Warning: \1", text, flags=re.MULTILINE)
    return text


# To demo capabilities
if __name__ == "__main__":
    text = """
<h1>Header 1</h1>
<h2>Header 2</h2>
<cmd>A command with a <parameter> or two</cmd>
<error>Error message</error>
<warning>Warning message</warning>
<success>Success message</success>
This is a link: <link>google.com</link>
<link>\x1b]8;;http://example.com\aThis is a masked link\x1b]8;;\a</link> (not supported in Mac Terminal)
- - -
<bold>Bold text</bold>
<soft>Soft text</soft>
<italic>Italic text</italic>
<underline>Underline text</underline>
<blink>Blinking text</blink>
<reverse>Reverse text</reverse>
Hidden text --><hidden>here</hidden><--
<strikethrough>Strikethrough text</strikethrough>
- - -
<black>Black text</black>
<red>Red text</red>
<green>Green text</green>
<yellow>Yellow text</yellow>
<blue>Blue text</blue>
<magenta>Magenta text</magenta>
<cyan>Cyan text</cyan>
<white>White text</white>
- - -
<bright_black>Bright black text</bright_black>
<bright_red>Bright red text</bright_red>
<bright_green>Bright green text</bright_green>
<bright_yellow>Bright yellow text</bright_yellow>
<bright_blue>Bright blue text</bright_blue>
<bright_magenta>Bright magenta text</bright_magenta>
<bright_cyan>Bright cyan text</bright_cyan>
<bright_white>Bright white text</bright_white>
- - -
<on_black>On black text</on_black>
<on_red>On red text</on_red>
<on_green>On green text</on_green>
<on_yellow>On yellow text</on_yellow>
<on_blue>On blue text</on_blue>
<on_magenta>On magenta text</on_magenta>
<on_cyan>On cyan text</on_cyan>
<on_white>On white text</on_white>
- - -
<on_bright_black>On bright black text</on_bright_black>
<on_bright_red>On bright red text</on_bright_red>
<on_bright_green>On bright green text</on_bright_green>
<on_bright_yellow>On bright yellow text</on_bright_yellow>
<on_bright_blue>On bright blue text</on_bright_blue>
<on_bright_magenta>On bright magenta text</on_bright_magenta>
<on_bright_cyan>On bright cyan text</on_bright_cyan>
<on_bright_white>On bright white text</on_bright_white>
"""
    print(__doc__)
    print_s(text, pad=2, tabs=1, edge=True, nowrap=True)
