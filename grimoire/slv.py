

SSOURCE = 'source'
STARGET = 'target'

class SLV(object):

    def __init__ ( self, coords, id, junction_list, type ):
        self.coords = coords
        self.id = id
        self.junction_list = junction_list
        if type != SSOURCE and type != STARGET: raise RuntimeError('Incorrect SLV type %s'%type)
        self.type = type

    def get_coordinates(self):
        return self.coords

    def get_junctions_list(self):
        return self.junction_list

    def is_Ssource(self):
        return bool(self.type==SSOURCE)

    def is_Starget(self):
        return bool(self.type == STARGET)
