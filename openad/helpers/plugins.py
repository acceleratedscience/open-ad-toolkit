def reorder_commands_by_category_index(plugin_commands: dict) -> list:
    """
    Reorder the commands by their index, per category.

    Parameters
    ----------
    plugin_commands : dict
        The plugin commands, unordered, eg:
        [
            { category: "Category A", command: "Foo", index: 1 },
            { category: "Category B", command: "Baz", index: 2 },
            { category: "Category A", command: "Bar", index: 0 },
            { category: "Category B", command: "Foo", index: 1 },
            { category: "Category A", command: "Baz", index: 2 },
            { category: "Category B", command: "Bar", index: 0 },
        ]

    Returns
    -------
    list
        The reordered plugin commands, eg:
        [
            { category: "Category A", command: "Bar", index: 0 },
            { category: "Category A", command: "Foo", index: 1 },
            { category: "Category A", command: "Baz", index: 2 },
            { category: "Category B", command: "Bar", index: 0 },
            { category: "Category B", command: "Foo", index: 1 },
            { category: "Category B", command: "Baz", index: 2 },
    """

    # Organize commands by category
    plugin_commands_organized = {}
    for Cmd in plugin_commands:
        if Cmd.category not in plugin_commands_organized:
            plugin_commands_organized[Cmd.category] = []
        plugin_commands_organized[Cmd.category].append(Cmd)

    # Reorder the commands by their index, per category
    plugin_commands_output = []
    for cat_cmds in plugin_commands_organized.values():
        cat_cmds.sort(key=lambda x: x.index)
        for Cmd in cat_cmds:
            plugin_commands_output.append(Cmd)

    return plugin_commands_output
