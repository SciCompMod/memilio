import os
from os.path import exists
from clang.cindex import *
from subprocess import check_output
import clang.cindex
import copy

class Scanner:
    def __init__(self, conf):
        self.conf = conf
        self.list_ast = []

        Config.set_library_file("/home/betz_mx/memilio/pycode/virtualenv/lib/python3.8/site-packages/libclang-14.0.1-py3.8-linux-x86_64.egg/clang/native/libclang.so")
        self.create_ast()
        #output_cursor_and_children(list_ast[1].cursor)

    def create_ast(self):
        """
        Creates an Ast for every .cpp/.h-file in the folder defined in the ScannerConf. 
        Saves them in list_ast.
        """
        # look if folder exists
        project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/models/"
        folder_path = project_path + self.conf.folder

        idx = Index.create()
        for _, _, files in os.walk(folder_path):
            for file_name in files:
                if file_name.endswith((".cpp", ".h")): 
                    ast = idx.parse(folder_path + "/" + file_name, args=['-x', 'c++'], options=0)
                    self.list_ast.append(ast)

    def extract_results(self, model):
        """
        Extracts the information of the asts and saves them in the data class model.
        Iterates over list of list_ast and calls find_node to visit all nodes of ast.
        """
        if type(self.list_ast) is not list:
            self.output_cursor(self.list_ast.cursor, 1)
            self.find_node(self.list_ast.cursor, model)
            return
        for ast in self.list_ast:
            self.output_cursor(ast.cursor, 1)
            self.find_node(ast.cursor, model)


    def find_node(self, node, model, namespace = ""):
        """
        Recursively walks over every node of an ast. Saves the namespace the node is in.
        Calls check_node_kind for extracting information out of the nodes.
        """
        if node.kind == CursorKind.NAMESPACE:
            namespace = (namespace + node.spelling + "::")
        else:
            self.check_node_kind(node, model, namespace)
        
        for n in node.get_children(): 
            self.find_node(n, model, namespace)
    
    def check_node_kind(self, node, model, namespace):
        """
        Checks the kind of node and calls the appropriate method for the given node kind.
        """
        kind = node.kind
        if node.is_scoped_enum(): # alternative kind == CursorKind.ENUM_DECL and node.spelling != "":
            model.enum_dict[namespace + node.spelling] = []
        elif kind == CursorKind.ENUM_CONSTANT_DECL and namespace + node.semantic_parent.spelling in model.enum_dict.keys():
            key = namespace + node.semantic_parent.spelling
            model.enum_dict[key].append(key + "::" + node.spelling)
        elif kind == CursorKind.CLASS_DECL:
            self.check_class(node, model, namespace)
        elif kind == CursorKind.CONSTRUCTOR:
            self.check_constructor(node, model, namespace)


    def check_class(self, node, model, namespace):
        """
        Inspect the nodes of kind CLASS_DECL and writes needed informations into model.
        Informations: model.name, model.population_groups
        """
        base_classes = self.get_base_class_string(node)
        if "CompartmentalModel" in base_classes:
            model.name = namespace + node.spelling
            for base_class in base_classes:
                if base_class.startswith("CompartmentalModel:Populations:"):
                    model.population_groups.append(base_class.lstrip("CompartmentalModel:Populations:"))
    
    def check_constructor(self, node, model, namespace):
        """
        Inspect the nodes of kind CONSTRUCTOR and writes needed informations into model.
        Informations: model.init
        """
        if (namespace + node.spelling) == model.name:
            init = {"type" : [], "name" : []}
            for arg in node.get_arguments():
                tokens = []
                for token in arg.get_tokens():
                    tokens.append(token.spelling)
                init["type"].append(" ".join(tokens[:-1]))
                init["name"].append(tokens[-1])
            model.init.append(init)

    @staticmethod
    def get_base_class_string(node):
        """
        Input: node of kind CLASS_DECL
        Returns the base class as a list of strings.
        
        Example: 
            Model : Base<Template> {} 
            returns ["Base", "Base:Template"]
        """

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

    def output_ast(self, index):
        """
        Outputs the ast at index in self.list_ast
        """
        self.output_cursor_and_children(self.list_ast[index].cursor)

    @classmethod
    def indent(cls, level):
        """ 
        Indentation string for pretty-printing
        """ 
        return '  '*level

    @classmethod
    def output_cursor(cls, cursor, level):
        """ 
        Low level cursor output
        """
        spelling = ''
        displayname = ''

        if cursor.spelling:
            spelling = cursor.spelling
        if cursor.displayname:
            displayname = cursor.displayname
        kind = cursor.kind

        print(cls.indent(level) + spelling, '<' + str(kind) + '>')
        print(cls.indent(level+1) + '"'  + displayname + '"')

    @classmethod
    def output_cursor_and_children(cls, cursor, level=0):
        """ 
        Output this cursor and its children with minimal formatting.
        """
        cls.output_cursor(cursor, level)
        if cursor.kind.is_reference():
            print(cls.indent(level) + 'reference to:')
            cls.output_cursor(cursor.referenced, level+1)

        # Recurse for children of this cursor
        has_children = False
        for c in cursor.get_children():
            if not has_children:
                print(cls.indent(level) + '{')
                has_children = True
            cls.output_cursor_and_children(c, level+1)

        if has_children:
            print(cls.indent(level) + '}')