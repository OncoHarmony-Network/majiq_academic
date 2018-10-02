import configparser
from os.path import join

from voila import constants


class Config:
    def __init__(self):
        # todo: this should be a singleton class
        self.files = None
        self.voila_files = None
        self.splice_graph_file = None
        self.read()

    @property
    def voila_file(self):
        return self.voila_files[0]

    @staticmethod
    def write(args):
        print('hey')
        config = configparser.ConfigParser()
        files = 'FILES'
        config.add_section(files)
        config.set(files, 'voila', '\n'.join(args.voila_files))
        config.set(files, 'splice_graph', args.splice_graph)
        with open(join(constants.CONFIG_FILE), 'w') as configfile:
            config.write(configfile)

    def read(self):
        config = configparser.ConfigParser()
        config.read(constants.CONFIG_FILE)
        self.files = config['FILES']
        self.voila_files = self.files['voila'].split('\n')
        self.splice_graph_file = self.files['splice_graph']
