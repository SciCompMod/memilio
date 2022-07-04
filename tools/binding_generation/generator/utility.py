from sympy.codegen.ast import Assignment, aug_assign

class Variable:
    def __init__(self, n, s):
        self.name = n
        self.symbol = s

class Parameter:
    def __init__(self, n, c, d, s):
        self.name = n
        self.c_type = c
        self.default = d
        self.symbol = s

def join(string_list, seperator=""):
    ret = string_list[0]
    for s in string_list[1:]:
        ret += seperator + s
    return ret

def assign(lhs, rhs, op="="):
    if op in ["", "=", ":="]:
        return Assignment(lhs, rhs)
    else:
        return aug_assign(lhs, op[0], rhs)
