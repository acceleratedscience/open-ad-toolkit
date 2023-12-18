# import pprint
# pp = pprint.PrettyPrinter(indent=4, width=50)


def hello_world(inputs: dict, cmd_pointer):
    reset = "\x1b[0m"
    red = "\x1b[31m"
    yellow = "\x1b[33m"
    green = "\x1b[32m"
    blue = "\x1b[34m"
    cyan = "\x1b[36m"
    magenta = "\x1b[35m"
    rainbow = [red, yellow, green, blue, cyan, magenta]
    inset = "   "

    # Hello, world rainbow
    print(f"\n{inset}", end="")
    for i, char in enumerate("hello, world"):
        color = rainbow[i % len(rainbow)]
        print(color + char, end="")
    print("\n")

    # inputs
    print(f"{inset}{yellow}inputs:{reset}", inputs, end="\n\n")

    # cmd_pointer
    print(f"{inset}{yellow}toolkit_current.toolkit_name:{reset}", cmd_pointer.toolkit_current.toolkit_name)
    print(f"{inset}{yellow}notebook_mode:{reset}", cmd_pointer.notebook_mode)
    print(f"{inset}{yellow}api_mode:{reset}", cmd_pointer.api_mode)
    print(f"{inset}{yellow}home_dir:{reset}", cmd_pointer.home_dir)
    print(f"{inset}{yellow}repo_dir:{reset}", cmd_pointer.repo_dir)
    print(f"{inset}{yellow}toolkit_dir:{reset}", cmd_pointer.toolkit_dir)
    print(f"{inset}{yellow}toolkit_current.methods:{reset}", cmd_pointer.toolkit_current.methods)
    print("")
