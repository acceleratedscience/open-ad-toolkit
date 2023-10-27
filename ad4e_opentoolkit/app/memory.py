import ad4e_opentoolkit.helpers.output as output

class Memory:
    def __init__(self):

        # Working memory
        # - - -
        # Here we store data which can be manipulated by follow-up commands.
        # If the next command is not one of the approved follow-up commands,
        # the memory is wiped to avoid performance impact.
        # - - -
        # Examples:
        # >> display data 'table1.csv' # Data is stored in short term memory.
        # >> result open # Data is used to spin up UI.
        self._storage = None

        # By default, the memory is wiped after each command,
        # unless the command is one of the approved follow-up commands,
        # in which case the memory is preserved via this dictionary.
        self._preserve = False

        # import random
        # self.test = random.randint(0,10)
    
    def __call__(self):
        print('call')
    

    def store(self, data):
        self._storage = data
        self._preserve = True
        # print(111, self.test)
    
    def get(self):
        # print(222, self.test)
        return self._storage
    
    def wipe(self):
        self._storage = None
    
    

    # This is run at the beginning of any follow up function
    # to check and preserve the memory content.
    def hold(self):
        # Preserve memory for further follow-ups.
        self._preserve = True

        # Abort if memory is empty.
        if self.store is None:
            output.output_error(output.msg('no_data_memory'), self)
    
    def holding(self):
        return self._preserve

    def release(self):
        self._preserve = False

memory = Memory()
