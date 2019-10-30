"""
Couple object: married pair of individuals and IBDs that one or both members is
expected to possess.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 8/13/19
"""
# DONE

class Couple:
    """
    Couple represents a married pair of Individuals and contains a list of
    IBDs that one or both member(s) of the couple is expected to possess.
    """

    def __init__(self, father, mother):
        self.p = father
        self.m = mother
        self.sex = -1 # changing to int (make hap unknown in descendants)
        self.id = self.p.id + "&" + self.m.id
        self._IBDs = {"00": set(), # unknown chromosome, uncertain
                    "01": set()    # unknown chromosome, certain
                }

    def __str__(self):
        return "p: %s, m: %s" %(self.p.id, self.m.id)

    def __eq__(self, other):
        return other != None and self.p.id == other.p.id and \
            self.m.id == other.m.id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.id)
