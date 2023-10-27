import ad4e_opentoolkit.helpers.output as output


class Memory:
    """
    Working memory
    --------------

    Description:
    The working memory is used to store data between commands.

    When data is stored in th memory, it can be used by follow-up commands.
    Follow up commands should always look like: `result <action>`.
    It is currently implemented by output_table() which we use to display data.
    The available follow-up commands like `result edit` are listed under the table.
    Please follow this template when adding new follow-up commands.

    Implementation:

        from ad4e_opentoolkit.app.memory import memory

        def main_command():
            # Store data in memory
            memory.store(table)

        def follow_up_command():
            # Don't erase the memory so it can be
            # reused by more follow-up commands.
            memory.preserve()

            # Retrieve data from memory
            const table = memory.get()

    Notes:
        - The memory is wiped after every command, unless memory.preserve() is called.
        - There are no restrictions as to what kind of data can be stored.

    In action:

        display data 'sample.csv'

        result open
    """

    def __init__(self):
        self._storage = None  # Where we store the data.
        self._preserve = False  # Inhibits memory wipe.

    #################
    # GETTERS
    #################

    def get(self):
        return self._storage

    #################
    # SETTERS
    #################

    # Preserve the memory for further follow up commands.
    # This prevents the memory from being wiped when the
    # command is executed. It should be placed at the
    # beginning of your command function.
    def preserve(self):
        self._preserve = True

        # Abort if memory is empty.
        if self.store is None:
            output.output_error(output.msg("no_data_memory"), self)

    def store(self, data):
        self._storage = data
        self._preserve = True

    def release(self):
        self._preserve = False

    def wipe(self):
        self._storage = None

    #################
    # FUNCTIONS
    #################

    def before_command(self):
        if memory._preserve:
            memory.release()  # Wipe with next command
        else:
            memory.wipe()  # Wipe now


memory = Memory()
