"""
Parse XML style tags for easy styling of text output.

Styles:
<h1>I'm a header</h1>
<cmd>ibmad ds foobar</cmd>

Colors:
<red>I'm red</red>

Supported colors names:
https://click.palletsprojects.com/en/8.1.x/api/#click.style
  • black (might be a gray)
  • red
  • green
  • yellow (might be an orange)
  • blue
  • magenta
  • cyan
  • white (might be light gray)
  • bright_red
  • bright_green
  • bright_yellow
  • bright_blue
  • bright_magenta
  • bright_cyan
  • bright_white

  • reset (reset the color code only)

"""

import re
import textwrap
from collections import deque

# We only look for allowed tags, this way we can still
# have xml type tags in the text without them being parsed.
# This is required because we use <these_tags> to denote
# command line arguments.
allowed_tags = [
    # Custom tags
    'h1',  # How we display primary headers: white with yellow underline
    'h2',  # How we display secondary headers: yellow
    'cmd',  # How we display commands: cyan
    'error',  # How we display errors: red
    'warning',  # How we display warnings: yellow
    'success',  # How we display success: green
    'link',  # How we display links: magenta (for now)
    'soft',  # Soft text: bright_black, effectively gray

    # Style tags
    'underline',
    'bold',

    # Color tags
    'reset',
    'black',
    'red',
    'green',
    'yellow',
    'blue',
    'magenta',
    'cyan',
    'white',
    'bright_black',
    'bright_red',
    'bright_green',
    'bright_yellow',
    'bright_blue',
    'bright_magenta',
    'bright_cyan',
    'bright_white',
    'on_black',
    'on_red',
    'on_green',
    'on_yellow',
    'on_blue',
    'on_magenta',
    'on_cyan',
    'on_white'
]


def style(text: str, pad: int = 0, pad_top: int = 0, xml_tags=True, pad_btm: int = 0, tabs=0, edge: bool | str = False, width=80, nowrap=False):
    """
    Parse xml tags and return styled output.

    Parameters:
    • pad (int): Number of line breaks before and after.
    • pad_top (int): Number of line breaks before.
    • pad_btm (int): Number of line breaks after.
    • tabs (int): Number of tabs to indent.
    • edge (bool|str): Line left edge with | Can be color or True for grey.
    • width (int): Width of paragraph in characters.
    • nowrap (bool): Do not wrap lines.

    Returns:
    Styled text ready to be printed to console.
    """

    if text is None:
        return ''

    # Enforce string type.
    text = str(text)

    # Strip any line breaks at start/end.
    text = _trim(text)

    # Add dotted lines under h1 tags.
    text = _add_header_lines(text)

    # Replace line breaks so all text is parsed on one line.
    text = text.replace('\n', '---LINEBREAK1---')

    # Make links clickable.
    text = _parse_links(text)

    # Style xml command parameters.
    text = style_cmd_params(text)

    # Style xml style tags.
    if xml_tags == True:
        text = parse_tags(text)

    # Restore line breaks.
    text = text.replace('---LINEBREAK1---', '\n')

    # Add tabs to the start of each line.
    if tabs or edge:
        tabs = tabs if tabs else 1
        text = add_tabs(text, tabs)

    # Wrap lines based on paragraph width in characters.
    if not nowrap:
        text = wrap_text(text, width=width)

    # Add edge.
    if edge:
        text = _edge(text, edge)

    # Add top and bottom padding
    if pad:
        padding = '\n' * pad
        text = padding + text + padding
    else:
        if pad_top:
            padding = '\n' * pad_top
            text = padding + text
        if pad_btm:
            padding = '\n' * pad_btm
            text = text + padding

    return text


def print_s(text: str, ephemeral=False, **kwargs):
    """
    Print styled text.

    Parameters:
        ephemeral (bool): Print next line on the same line. (Currently not used)
    """

    text = str(text)
    end = '\r' if ephemeral else '\n'
    print(style(text, **kwargs), end=end)


def add_tabs(text, tabs):
    """ Add left padding to a text block. """

    text = str(text)
    spacing = '    ' * tabs
    return spacing + ('\n' + spacing).join(text.split('\n'))


def wrap_text(text: str, width=80):
    """ Wrap lines to fit within the paragraph width. """

    text = str(text)
    lines = text.splitlines()
    wrapped_lines = []

    for line in lines:
        stripped_line = line.lstrip()
        indent = line[:len(line) - len(stripped_line)]
        wrapped = textwrap.fill(stripped_line, width=width, subsequent_indent=indent)
        wrapped_lines.append(indent + wrapped)

    return '\n'.join(wrapped_lines)


def style_cmd_params(text: str):
    """
    Mark non-style related tags in bright_black.

    Some context:
    We use style tags, eg. <red>red text</red>, to style text.
    These style tags get replaced with style characters in parse_tags.
    We also display command parameters in the same format, eg. run <run_name>
    These parameters we display in grey, making it easier to read the commands.

    Because they always appear inside the <cmd> tag, we first need to close
    the cmd tag, then open the bright_black tag, then open the cmd tag again.

    So in practice:
    <cmd>run <run_name></cmd> --> <cmd>run </cmd><bright_black><run_name></bright_black><cmd>
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
            pattern = fr'<cmd>.*?(<(?!{"|".join(allowed_tags)})[^</>]+>).*</cmd>'
            param = re.findall(pattern, line_copy)
            if param:
                param = param[0]
                params.append(param)
                line_copy = line_copy.replace(param, '')
                line = line.replace(
                    param, f'<bright_black>{param}</bright_black>')
            else:
                searching = False

        result.append(line)

    return '\n'.join(result)


def parse_tags(text: str):
    """ Parse xml tags and return styled output. """

    # Enforce string type.
    text = str(text)

    # Remove <pre> tags which are used to wrap ASCII headers (see splash.py)
    text = re.sub(r"<pre[^<]*?>(.*?)</pre>", r'\1', text)  # flags=re.MULTILINE

    # Replace open-close tags <red>lorem ipsum</red>.
    pattern = fr"(.*?)<({'|'.join(allowed_tags)})>(.*?)</\2>"
    text = re.sub(pattern, lambda match: _replace(match, pattern), text)

    # Replace singular tags <red/> --> everything red from here on (not used).
    pattern_singular = r"<({})/>".format('|'.join(allowed_tags))
    text = re.sub(pattern_singular, _replace_singular, text)
    return text


def strip_tags(text: str):
    """ Recursively remove all XML tags. """

    def _strip(match: object, pattern):
        """ Replace regex matches with appropriate styling. """
        inner_text = match.group(2)

        # Strip any nested tags.
        if re.findall(pattern, inner_text):
            inner_text = re.sub(pattern, lambda match: _strip(match, pattern), inner_text)

        return inner_text

    # Enforce string type.
    text = str(text)

    # Replace line breaks so all text is parsed on one line.
    text = text.replace('\n', '---LINEBREAK2---')

    # Expand error and success tags.
    text = _expand_error_success_tags(text, False)

    # Strip tags.
    pattern = fr"<({'|'.join(allowed_tags)})>(.*?)</\1>"
    text = re.sub(pattern, lambda match: _strip(match, pattern), text)

    # Restore line breaks.
    text = text.replace('---LINEBREAK2---', '\n')

    return text


def tags_to_markdown(text: str):
    if text is None:
        return ''

    # Enforce string type.
    text = str(text)

    # Replace leading spaces with non-breaking spaces.
    text = re.sub(r'^( *)', '', text, flags=re.MULTILINE)

    # Replace line breaks so all text is parsed on one line.
    # Because html breaks (<br>) don't play well with headings,
    # and end of line characters don't play well with `code`
    # blocks, we have to do some trickery here.
    text = re.sub(r'(</h[123]>)(\n+)', lambda match: match.group(1) + len(match.group(2)) * '---LINEBREAKSOFT---', text)
    text = re.sub(r'(\n+)(<h[123]>)', lambda match: len(match.group(1)) * '---LINEBREAKSOFT---' + match.group(2), text)
    text = text.replace('\n', '---LINEBREAK3---')

    # Expand error and success tags.
    text = _expand_error_success_tags(text, True)

    # Replace tags
    # We only replace <soft> and <underline> tags
    # when they don't appear inside <cmd> tags.
    text = re.sub(r'(?<!\<cmd\>)<soft>([^<]*)</soft>(?!\</cmd\>)', r'<span style="color: #ccc">\1</span>', text)
    text = re.sub(r'(?<!\<cmd\>)<underline>([^<]*)<\/underline>(?!\</cmd\>)', r'<span style="text-decoration: underline">\1</span>', text)
    text = re.sub(r'<h1>(.*?)<\/h1>', r'## \1', text)
    text = re.sub(r'<h2>(.*?)<\/h2>', r'### \1', text)
    text = re.sub(r'<link>(.*?)<\/link>', r'<a target="_blank" href="\1">\1</a>', text)
    text = re.sub(r'<bold>(.*?)<\/bold>', r'**\1**', text)
    text = re.sub(r'<cmd>(.*?)<\/cmd>', r'`\1`', text)
    text = re.sub(r'<on_red>(.*?)<\/on_red>', r'<span style="background: #d00; color: #fff">\1</span>', text)
    text = re.sub(r'<on_green>(.*?)<\/on_green>', r'<span style="background: #0d0; color: #fff">\1</span>', text)

    # Escape quotes.
    text = text.replace("'", "\'")

    # Replace all other tags
    text = strip_tags(text)

    # Restore line breaks.
    text = text.replace('---LINEBREAKSOFT---', '\n')
    text = text.replace('---LINEBREAK3---', '<br>')

    # TRASH
    # Sometimes we need to use a hard break – see adccl_intro.
    # text = text.replace('\n•', '<br>•')  # For bullet points
    # text = text.replace('\n`', '<br>`')  # When listing commands
    # text = text.replace('\n<', '<br><')  # Consecutive HTML elements
    # text = text.replace('>\n', '><br>')  # HTML element followed by regular text

    return text

#
#
#
#


# Color codes
_colors = {
    'underline': '\u001b[4m',
    'bold': '\u001b[1m',
    'reset': '\u001b[0m',

    # Foreground colors - Regular
    'black': '\u001b[30m',
    'red': '\u001b[31m',
    'green': '\u001b[32m',
    'yellow': '\u001b[33m',
    'blue': '\u001b[34m',
    'magenta': '\u001b[35m',
    'cyan': '\u001b[36m',
    'white': '\u001b[37m',

    # Forground colors - Bright (non-standard)
    'bright_black': '\u001b[90m',
    'bright_red': '\u001b[91m',
    'bright_green': '\u001b[92m',
    'bright_yellow': '\u001b[93m',
    'bright_blue': '\u001b[94m',
    'bright_magenta': '\u001b[95m',
    'bright_cyan': '\u001b[96m',
    'bright_white': '\u001b[97m',

    # Background colors - Regular
    'on_black': '\u001b[0;40m',
    'on_red': '\u001b[0;41m',
    'on_green': '\u001b[0;42m',
    'on_yellow': '\u001b[0;43m',
    'on_blue': '\u001b[0;44m',
    'on_magenta': '\u001b[0;45m',
    'on_cyan': '\u001b[0;46m',
    'on_white': '\u001b[0;47m',

    # Foreground colors - Bright (non-standard)
    'on_bright_black': '\u001b[100m',
    'on_bright_red': '\u001b[101m',
    'on_bright_green': '\u001b[102m',
    'on_bright_yellow': '\u001b[103m',
    'on_bright_blue': '\u001b[104m',
    'on_bright_magenta': '\u001b[105m',
    'on_bright_cyan': '\u001b[106m',
    'on_bright_white': '\u001b[107m',

    # Unused but useful as a reference
    # \u001b[7m - Reversed
}


def _trim(text: str):
    """ Remove line breaks at front/end caused by putting quotes on separate line. """
    lines = deque(text.split('\n'))

    if lines[0].strip() == '':
        lines.popleft()
    if lines[len(lines) - 1].strip() == '':
        lines.pop()

    return '\n'.join(lines)


def _expand_error_success_tags(text, html=False):
    """ Turn error/success tags into text, for non-styles output. """
    if html:
        # text = re.sub(r'^<error>(.*)<\/error>$', r'<span style="color: red">Error: \1</span>', text)
        # text = re.sub(r'^<success>(.*)<\/success>$', r'<span style="color: green">Success: \1</span>', text)
        # text = re.sub(r'^<warning>(.*)<\/warning>$', r'<span style="color: orange">Warning: \1</span>', text)
        text = re.sub(r'^<error>(.*)<\/error>', r'<span style="color: red">\1</span>', text, flags=re.MULTILINE)
        text = re.sub(r'^<success>(.*)<\/success>', r'<span style="color: green">\1</span>', text, flags=re.MULTILINE)
        text = re.sub(r'^<warning>(.*)<\/warning>', r'<span style="color: orange">\1</span>', text, flags=re.MULTILINE)
    else:
        text = re.sub(r'^<error>(.*)<\/error>', r'Error: \1', text, flags=re.MULTILINE)
        text = re.sub(r'^<success>(.*)<\/success>', r'Success: \1', text, flags=re.MULTILINE)
        text = re.sub(r'^<warning>(.*)<\/warning>', r'Warning: \1', text, flags=re.MULTILINE)
    return text


def _add_header_lines(text):
    """ Add dotted lines under h1 tags."""
    pattern = fr"(.*?)<h1>(.*?)</h1>"

    def _modify(match: object):
        original_text = match.group(0)
        space_before = match.group(1)
        inner_text = match.group(2)
        lines = space_before + '<yellow>' + ('-' * len(strip_tags(inner_text))) + '</yellow>'
        return original_text + '\n' + lines

    return re.sub(pattern, _modify, text)


def _parse_links(text):
    """ Parse html links and return clickable output. """
    lines = text.splitlines()
    result = []
    for line in lines:

        # Detect <a href="x.yz">click here</a>
        pattern = r'<a href="(.*)">(.*?)<\/a>'
        link_tuples = re.findall(pattern, line)  # [(url, text), ...]
        match_strings = re.search(pattern, line)  # ['<a href="xyz">click here</a>', ...]
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

    return '\n'.join(result)


def _replace(match: object, pattern, parent_color='reset'):
    """ Replace regex matches with appropriate styling. """
    text_before = match.group(1)
    tag = match.group(2)
    inner_text = match.group(3)
    color = tag_to_color(tag)

    parent_color = tag_to_color(parent_color)
    color_code_open = _colors[color]
    color_code_close = _colors[parent_color]
    # We need reset_code before color_code_clode
    # for underline/bold to work.
    reset_code = _colors['reset']

    # Replace any nested tags.
    if re.findall(pattern, inner_text):
        inner_text = re.sub(pattern, lambda match: _replace(match, pattern, color), inner_text)

    return f"{text_before}{color_code_open}{inner_text}{reset_code}{color_code_close}"


def _replace_singular(match: object):
    tag = match.group(1)
    color = tag_to_color(tag)
    return _colors[color]


def tag_to_color(tag):
    if tag == 'h1':
        color = 'reset'
    elif tag == 'h2':
        color = 'yellow'
    elif tag == 'cmd':
        color = 'cyan'
    elif tag == 'error':
        color = 'red'
    elif tag == 'warning':
        color = 'yellow'
    elif tag == 'success':
        color = 'green'
    elif tag == 'link':
        color = 'magenta'
        # color = 'underline'
    elif tag == 'soft':
        color = 'bright_black'
    elif tag in allowed_tags:
        color = tag
    else:
        color = 'reset'
    return color


def _edge(text, edge):
    """ Add edge to left side of a text block. """
    edge_color = 'bright_black' if edge == True else edge
    edge_str = '|'
    lines = text.splitlines()
    placc = ''  # Previous line ANSI color character
    for i, line in enumerate(lines):
        edge_str = _colors[edge_color] + edge_str + placc
        placc = re.findall(r'\x1b\[\d+m', line)
        placc = placc[len(placc) - 1] if len(placc) > 0 else ''
        lines[i] = edge_str + line

    return '\n'.join(lines)


# To demo capabilities
if __name__ == "__main__":
    text = """
<h1>Hello World</h1>
<cmd>ibmad ds foobar</cmd>

<black>Black text</black>
<red>Red text</red>
<green>Green text</green>
<yellow>Yellow text</yellow>
<blue>Blue text</blue>
<magenta>Magenta text</magenta>
<cyan>Cyan text</cyan>
<white>White text</white>

<bright_black>Bright black text</bright_black>
<bright_red>Bright red text</bright_red>
<bright_green>Bright green text</bright_green>
<bright_yellow>Bright yellow text</bright_yellow>
<bright_blue>Bright blue text</bright_blue>
<bright_magenta>Bright magenta text</bright_magenta>
<bright_cyan>Bright cyan text</bright_cyan>
<bright_white>Bright white text</bright_white>

<on_black>On black text</on_black>
<on_red>On red text</on_red>
<on_green>On green text</on_green>
<on_yellow>On yellow text</on_yellow>
<on_blue>On blue text</on_blue>
<on_magenta>On magenta text</on_magenta>
<on_cyan>On cyan text</on_cyan>
<on_white>On white text</on_white>

<on_bright_black>On bright black text</on_bright_black>
<on_bright_red>On bright red text</on_bright_red>
<on_bright_green>On bright green text</on_bright_green>
<on_bright_yellow>On bright yellow text</on_bright_yellow>
<on_bright_blue>On bright blue text</on_bright_blue>
<on_bright_magenta>On bright magenta text</on_bright_magenta>
<on_bright_cyan>On bright cyan text</on_bright_cyan>
<on_bright_white>On bright white text</on_bright_white>



"""
    styled_text = style(text, pad=2, tabs=1, edge=True)
    print(styled_text)
    print(_colors['green'] + 'ABCDEFGH' + _colors['reset'])
    print(_colors['blue'] + 'ABCDEFGH' + _colors['reset'])
