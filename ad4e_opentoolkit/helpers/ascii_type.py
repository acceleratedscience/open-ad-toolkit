"""ASCII title generator."""

# autopep8: off

abc = [
    """
  ░░░░  
░░    ░░
░░░░░░░░
░░    ░░
░░    ░░
""",
    """
░░░░░░  
░░    ░░
░░░░░░  
░░    ░░
░░░░░░░░
""",
    """
  ░░░░░░
░░      
░░      
░░      
  ░░░░░░
""",
    """
░░░░░░  
░░    ░░
░░    ░░
░░    ░░
░░░░░░  
""",
    """
░░░░░░░░
░░      
░░░░░░  
░░      
░░░░░░░░
""",
    """
░░░░░░░░
░░      
░░░░░░  
░░      
░░      
""",
    """
  ░░░░░░
░░      
░░  ░░░░
░░    ░░
  ░░░░░░
""",
    """
░░    ░░
░░    ░░
░░░░░░░░
░░    ░░
░░    ░░
""",
    """
░░░░░░
  ░░  
  ░░  
  ░░  
░░░░░░
""",
    """
  ░░░░░░
      ░░
      ░░
░░    ░░
  ░░░░  
""",
    """
░░    ░░
░░  ░░  
░░░░    
░░  ░░  
░░    ░░
""",
    """
░░      
░░      
░░      
░░      
░░░░░░░░
""",
    """
░░      ░░
░░░░  ░░░░
░░  ░░  ░░
░░      ░░
░░      ░░
""",
    """
░░      ░░
░░░░    ░░
░░  ░░  ░░
░░    ░░░░
░░      ░░
""",
    """
  ░░░░░░  
░░      ░░
░░      ░░
░░      ░░
  ░░░░░░  
""",
    """
░░░░░░  
░░    ░░
░░    ░░
░░░░░░  
░░      
""",
    """
  ░░░░░░
░░    ░░
░░    ░░
  ░░░░░░
      ░░
""",
    """
░░░░░░  
░░    ░░
░░    ░░
░░░░░░  
░░    ░░
""",
    """
  ░░░░░░
░░      
  ░░░░  
      ░░
░░░░░░  
""",
    """
░░░░░░░░░░
    ░░    
    ░░    
    ░░    
    ░░    
""",
    """
░░    ░░
░░    ░░
░░    ░░
░░    ░░
  ░░░░  
""",
    """
░░      ░░
░░      ░░
░░      ░░
  ░░  ░░  
    ░░    
""",
    """
░░      ░░
░░  ░░  ░░
░░  ░░  ░░
  ░░  ░░  
  ░░  ░░  
""",
    """
░░      ░░
  ░░  ░░  
    ░░    
  ░░  ░░  
░░      ░░
""",
    """
░░      ░░
░░      ░░
  ░░  ░░  
    ░░    
    ░░    
""",
    """
░░░░░░░░
      ░░
  ░░░░  
░░      
░░░░░░░░
""",
    """
  ░░░░  
░░    ░░
░░░░░░░░
░░    ░░
  ░░░░  
""",
    """
  ░░░░  
    ░░  
    ░░  
    ░░  
░░░░░░░░
""",
    """
  ░░░░  
░░    ░░
    ░░  
  ░░    
░░░░░░░░
""",
    """
  ░░░░  
░░    ░░
    ░░░░
░░    ░░
  ░░░░  
""",
    """
░░    ░░
░░    ░░
  ░░░░░░
      ░░
      ░░
""",
    #     '''
    #       ░░
    #     ░░░░
    #   ░░  ░░
    # ░░░░░░░░
    #       ░░
    # ''',
    """
░░░░░░░░
░░      
░░░░░░  
      ░░
░░░░░░  
""",
    """
  ░░░░  
░░      
░░░░░░  
░░    ░░
  ░░░░  
""",
    """
░░░░░░░░
      ░░
    ░░  
  ░░    
  ░░    
""",
    """
  ░░░░  
░░    ░░
  ░░░░  
░░    ░░
  ░░░░  
""",
    """
  ░░░░  
░░    ░░
  ░░░░░░
      ░░
  ░░░░  
""",
    " ",
]

abc_reverse = [
    # A
    """
░░    ░░
  ░░░░  
        
  ░░░░  
  ░░░░  
""",
    # B
    """
      ░░
  ░░░░  
      ░░
  ░░░░  
        
""",
    # C
    """
░░      
  ░░░░░░
  ░░░░░░
  ░░░░░░
░░      
""",
    # D
    """
      ░░
  ░░░░  
  ░░░░  
  ░░░░  
      ░░
""",
    # E
    """
        
  ░░░░░░
      ░░
  ░░░░░░
        
""",
    # F
    """
        
  ░░░░░░
      ░░
  ░░░░░░
  ░░░░░░
""",
    # G
    """
░░      
  ░░░░░░
  ░░    
  ░░░░  
░░      
""",
    # H
    """
  ░░░░  
  ░░░░  
        
  ░░░░  
  ░░░░  
""",
    # I
    """
      
░░  ░░
░░  ░░
░░  ░░
      
""",
    # J
    """
░░      
░░░░░░  
░░░░░░  
  ░░░░  
░░    ░░
""",
    # K
    """
  ░░░░  
  ░░  ░░
    ░░░░
  ░░  ░░
  ░░░░  
""",
    # L
    """
  ░░░░░░
  ░░░░░░
  ░░░░░░
  ░░░░░░
        
""",
    # M
    """
  ░░░░░░  
    ░░    
  ░░  ░░  
  ░░░░░░  
  ░░░░░░  

""",
    # N
    """
  ░░░░░░  
    ░░░░  
  ░░  ░░  
  ░░░░    
  ░░░░░░  
""",
    # O
    """
░░      ░░
  ░░░░░░  
  ░░░░░░  
  ░░░░░░  
░░      ░░
""",
    # P
    """
      ░░
  ░░░░  
  ░░░░  
      ░░
  ░░░░░░
""",
    # Q
    """
░░      
  ░░░░  
  ░░░░  
░░      
░░░░░░  
""",
    # R
    """
      ░░
  ░░░░  
  ░░░░  
      ░░
  ░░░░  
""",
    # S
    """
░░      
  ░░░░░░
░░    ░░
░░░░░░  
      ░░
""",
    # T
    """
          
░░░░  ░░░░
░░░░  ░░░░
░░░░  ░░░░
░░░░  ░░░░
""",
    # U
    """
  ░░░░  
  ░░░░  
  ░░░░  
  ░░░░  
░░    ░░
""",
    # V
    """
  ░░░░░░  
  ░░░░░░  
  ░░░░░░  
░░  ░░  ░░
░░░░  ░░░░
""",
    # W
    """
  ░░░░░░  
  ░░  ░░  
  ░░  ░░  
░░  ░░  ░░
░░  ░░  ░░
""",
    # X
    """
  ░░░░░░  
░░  ░░  ░░
░░░░  ░░░░
░░  ░░  ░░
  ░░░░░░  
""",
    # Y
    """
  ░░░░░░  
  ░░░░░░  
░░  ░░  ░░
░░░░  ░░░░
░░░░  ░░░░
""",
    # Z
    """
        
░░░░░░  
░░    ░░
  ░░░░░░
        
""",
    # 0
    """
░░    ░░
  ░░░░  
        
  ░░░░  
░░    ░░
""",
    # 1
    """
░░    ░░
░░░░  ░░
░░░░  ░░
░░░░  ░░
        
""",
    # 2
    """
░░    ░░
  ░░░░  
░░░░  ░░
░░  ░░░░
        
""",
    # 3
    """
░░    ░░
  ░░░░  
░░░░    
  ░░░░  
░░    ░░
""",
    # 4
    """
  ░░░░  
  ░░░░  
░░      
░░░░░░  
░░░░░░  
""",
    # 5
    """
        
  ░░░░░░
      ░░
░░░░░░  
      ░░
""",
    # 6
    """
░░    ░░
  ░░░░░░
      ░░
  ░░░░  
░░    ░░
""",
    # 7
    """
        
░░░░░░  
░░░░  ░░
░░  ░░░░
░░  ░░░░
""",
    # 8
    """
░░    ░░
  ░░░░  
░░    ░░
  ░░░░  
░░    ░░
""",
    # 9
    """
░░    ░░
  ░░░░  
░░      
░░░░░░  
░░    ░░
""",
    # Space
    """
░░
░░
░░
░░
░░
""",
]

CHARS = "abcdefghijklmnopqrstuvwxyz0123456789 "
abc = {letter: abc[CHARS.index(letter)].strip("\n") for letter in CHARS}
abc_reverse = {letter: abc_reverse[CHARS.index(letter)].strip("\n") for letter in CHARS}


def ascii_type(text, reverse=False, char=None):
    """
    Convert a string to big ASCII text.
    Supports a-z0-9 and space, case-insensitive.

    Parameters
    ----------
    text : str
        The text to display.
    reverse : bool, optional
        Reverses the foreground and background.
    char : str|int, optional
        The character used to paint.
    """

    # Parse text ans set style.
    chars = ["░", "█"]
    the_abc = abc_reverse if reverse else abc
    ascii_letters = [the_abc[letter] for letter in text.lower()]
    if char is not None:
        if isinstance(char, int):
            char = chars[char]
        ascii_letters = list(map(lambda x: x.replace("░", char), ascii_letters))
    else:
        char = chars[0]
    blank = char * 2 if reverse else "  "

    # Measure total line length.
    line_len = 4 if reverse else 0
    for i, letter in enumerate(ascii_letters):
        first_line = letter.splitlines()[0]
        line_len += len(first_line)
        if i < len(ascii_letters) - 1:
            line_len += 2

    output = ""

    # Top padding for reverse text.
    if reverse:
        output += (blank * int(line_len / 2)) + "\n"

    for line in range(5):
        # Left padding for reverse text.
        if reverse:
            output += blank

        # Compile line.
        for i, letter in enumerate(ascii_letters):
            # Space
            if letter == " ":
                output += blank

            # Letter
            else:
                output += letter.split("\n")[line]
                # Space between letters.
                if i < len(ascii_letters) - 1:
                    output += blank

        # Right padding for reverse text.
        if reverse:
            output += blank

        output += "\n"

    # Bottom padding for reverse text.
    if reverse:
        output += (blank * int(line_len / 2)) + "\n"

    text = "".join(output)
    return text
