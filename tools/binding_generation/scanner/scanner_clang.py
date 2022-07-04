import os
from os.path import exists
from clang.cindex import *
import clang.cindex
import copy

def get_infos_clang(model, conf):
    
    Config.set_library_file("/home/betz_mx/memilio/pycode/virtualenv/lib/python3.8/site-packages/libclang-14.0.1-py3.8-linux-x86_64.egg/clang/native/libclang.so")

    # look if folder exists
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/models/"
    folder_path = project_path + conf.folder
    
    
    # scanner
    list_tu = get_ast(folder_path)
    extract_results(list_tu, model)
    #output_cursor_and_children(list_tu[1].cursor)
    print(model)



def extract_results(list_tu, model):
    if type(list_tu) is not list:
        output_cursor(list_tu.cursor, 1)
        find_node(list_tu.cursor, model)
        return
    for ast in list_tu:
        output_cursor(ast.cursor, 1)
        find_node(ast.cursor, model)

def find_node(node, model, namespace = ""):
    if node.kind == CursorKind.NAMESPACE:
        namespace = (namespace + node.spelling + "::")
    else:
        check_node(node, model, namespace)
    
    for n in node.get_children(): 
        find_node(n, model, namespace)
    

def check_node(node, model, namespace):

    """
    Interesting
    """
    kind = node.kind
    # save enums
    if node.is_scoped_enum(): # alternative kind == CursorKind.ENUM_DECL and node.spelling != "":
        model.enum_dict[namespace + node.spelling] = []
    elif kind == CursorKind.ENUM_CONSTANT_DECL and namespace + node.semantic_parent.spelling in model.enum_dict.keys():
        key = namespace + node.semantic_parent.spelling
        model.enum_dict[key].append(key + "::" + node.spelling)
    elif kind == CursorKind.CLASS_DECL:
        base_classes = get_base_class(node)
        if "CompartmentalModel" in base_classes:
            model.name = namespace + node.spelling
            for base_class in base_classes:
                if base_class.startswith("CompartmentalModel:Populations:"):
                    model.population_groups.append(base_class.lstrip("CompartmentalModel:Populations:"))
    elif kind == CursorKind.CONSTRUCTOR:
        if (namespace + node.spelling) == model.name:
            init = {"type" : [], "name" : []}
            for arg in node.get_arguments():
                tokens = []
                for token in arg.get_tokens():
                    tokens.append(token.spelling)
                init["type"].append(" ".join(tokens[:-1]))
                init["name"].append(tokens[-1])
            model.init.append(init)

def get_base_class(node):

    if node.kind != CursorKind.CLASS_DECL:
        return None
    prefix = ""
    is_base = False
    result = []

    for token in node.get_tokens():
        
        if token.spelling == "{" or token.spelling == ";":
                return result 
        if not is_base:
            if token.spelling == ":":
                is_base = True
            continue

        if token.kind == TokenKind.IDENTIFIER:
            result.append((prefix + ":" + token.spelling).lstrip(":"))
        elif token.kind == TokenKind.PUNCTUATION:
            if token.spelling == "<":
                prefix = result[-1]
            elif token.spelling == ">":
                prefix = ":".join(prefix.split(":")[:-1])
    return result


def get_ast(folder_path):
    idx = Index.create()
    list_tu = []
    for _, _, files in os.walk(folder_path):
        for file_name in files:
            if file_name.endswith((".cpp", ".h")): 
                list_tu.append(idx.parse(folder_path + "/" + file_name, args=['-x', 'c++'], options=0))
    return list_tu


def indent(level):
    """ Indentation string for pretty-printing
    """ 
    return '  '*level

def output_cursor(cursor, level):
    """ Low level cursor output
    """
    spelling = ''
    displayname = ''

    if cursor.spelling:
        spelling = cursor.spelling
    if cursor.displayname:
        displayname = cursor.displayname
    kind = cursor.kind

    print(indent(level) + spelling, '<' + str(kind) + '>')
    print(indent(level+1) + '"'  + displayname + '"')

def output_cursor_and_children(cursor, level=0):
    """ Output this cursor and its children with minimal formatting.
    """
    output_cursor(cursor, level)
    if cursor.kind.is_reference():
        print(indent(level) + 'reference to:')
        #output_cursor(clang.cindex.Cursor_ref(cursor), level+1)

    # Recurse for children of this cursor
    has_children = False
    for c in cursor.get_children():
        if not has_children:
            print(indent(level) + '{')
            has_children = True
        output_cursor_and_children(c, level+1)

    if has_children:
        print(indent(level) + '}')