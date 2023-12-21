def test_clauses1(inputs: dict, cmd_pointer):
    reset = "\x1b[0m"
    yellow = "\x1b[33m"
    inset = "   "
    print(f"{inset}{yellow}inputs:{reset}", inputs, end="\n\n")

    if "test" in inputs:
        print('IN_CLAUSE clause detected --> "test"')

    if "foo" in inputs:
        print('USING clause detected --> "foo"')

    if "bar" in inputs:
        print('USING clause detected --> "bar"')

    if "baz" in inputs:
        print('USING clause detected --> "baz"')

    if "save_as" in inputs:
        print("SAVE_AS clause detected")

    if "estimate_only" in inputs:
        print("ESTIMATE_ONLY clause detected")

    if "return_as_data" in inputs:
        print("RETURN_AS_DATA clause detected")
    
    if "show_foo" in inputs:
        for option in inputs["show_foo"]:
            if option == "foo":


# def demo_clause2(inputs: dict, cmd_pointer):
#     reset = "\x1b[0m"
#     yellow = "\x1b[33m"
#     inset = "   "
#     print(f"{inset}{yellow}inputs:{reset}", inputs, end="\n\n")

#     if "save_as" in inputs:
#         print("SAVE_AS clause detected")

#     if "estimate_only" in inputs:
#         print("ESTIMATE_ONLY clause detected")

#     if "return_as_data" in inputs:
#         print("RETURN_AS_DATA clause detected")
