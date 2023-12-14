"Result Memory Library"
import openad.helpers.output as output


class Memory:
    """
    Working memory
    --------------

    Description:
    The working memory is used to store data between commands.

    When data is stored in the memory, it can be used by the follow-up commands.
    The follow up commands are structured as: `result <action>`.

    Currently data is stored by output_table() which we use to display data.
    The available follow-up commands like `result edit` are listed under the table.

    Implementation:

        def main_command(cmd_pointer):
            # Store data in memory (this happens via output_table() now)
            cmd_pointer.memory.store(table)

        def follow_up_command(cmd_pointer):
            # Don't erase the memory so it can be
            # reused by more follow-up commands.
            cmd_pointer.memory.preserve()

            # Retrieve data from memory
            const table = cmd_pointer.memory.get()

    Notes:
        - The memory is wiped after every command, unless memory.preserve() is called.
        - Follow-up commands expect table data, but we could store any type of data in the memory.

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
    # METHODS
    #################

    def before_command(self):
        if self._preserve:
            self.release()  # Wipe with next command
        else:
            self.wipe()  # Wipe now
