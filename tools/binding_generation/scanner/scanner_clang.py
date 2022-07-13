import os
from os.path import exists
from clang.cindex import *
from subprocess import check_output, call
import clang.cindex
import copy
from generator.Model import Model

class Scanner:
    def __init__(self, conf):

        Config.set_library_file(conf["scanner_input"]["libclang_library_path"])
        self.folder = conf["scanner_input"]["source_folder"]
        self.parser_information = conf["parser_information"]
        self.data_information = conf["data_information"]
        
        #call("")
        self.list_ast = []
        if conf["database"]["use_database"]:
            self.create_database_ast(conf["database"]["main_file"])
        else:
            self.create_ast()
        for dig in self.list_ast[0].diagnostics:
            print(dig)
        #output_cursor_and_children(list_ast[1].cursor)

    def create_ast(self):
        """
        Creates an Ast for every .cpp/.h-file in the folder defined in the ScannerConf. 
        Saves them in list_ast.
        """
        # look if folder exists
        project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/models/"
        folder_path = project_path + self.folder
        idx = Index.create()
        file_args = ['CMakeFiles/ode_seir.dir/model.cpp.o', "-c", "-Xclang", "-ast-dump", "-fsyntax-only"]
        ast = idx.parse(folder_path + "/" + "test.cpp", args=file_args)
        self.list_ast.append(ast)

    def create_ast2(self):
        """
        Creates an Ast for every .cpp/.h-file in the folder defined in the ScannerConf. 
        Saves them in list_ast.
        """
        # look if folder exists
        project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/models/"
        folder_path = project_path + self.folder
        idx = Index.create()
        file_args = ["clang++", 'CMakeFiles/ode_seir.dir/model.cpp.o', "-c", "-Xclang", "-ast-dump", "-fsyntax-only"]
        for _, _, files in os.walk(folder_path):
            for file_name in files:
                if file_name.endswith((".cpp")):#, ".h")): 
                    #ast = idx.parse(folder_path + "/" + file_name, args=['-x', 'c++'], options=0)
                    ast = idx.parse(folder_path + "/" + file_name, file_args)
                    self.list_ast.append(ast)

    def create_database_ast(self, file_name):
        """
        Creates an Ast for a main .cpp file with database. Requires an compile_commands.json.
        Saves them in list_ast.
        """
        # look if folder exists
        project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/models/"
        folder_path = project_path + self.folder
        idx = Index.create()
        
        # Step 1: load the compilation database
        compdb = CompilationDatabase.fromDirectory(check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/build")
        commands = compdb.getCompileCommands(folder_path + "/" + file_name)
        file_args = []
        for command in commands:
          for argument in command.arguments:
                if argument != '-DKF_IN_UE=1' and argument != '-DWITH_EDITOR=1':
                    # if '-IG:\\RedApp\\Plugins\\' in argument:
                    file_args.append(argument)

        file_args = file_args[0:-1]
        #file_args.append(r'-DDLLIMPORT=')
        #file_args.append(r'-DDLLEXPORT=')
        file_args.append(r'-Xclang')
        #file_args.append(r'-ast-dump')
        #file_args.append(r'-Wdeprecated')
        file_args.append(r'-fsyntax-only')
        #file_args.append(r'-Wunused-command-line-argument')
        ast = idx.parse(folder_path + "/" + "test.cpp", file_args)
        self.list_ast.append(ast)      
        
    def extract_results(self):
        """
        Extracts the information of the asts and saves them in the data class model.
        Iterates over list of list_ast and calls find_node to visit all nodes of ast.
        """
        model = Model()
        if type(self.list_ast) is not list:
            output_cursor(self.list_ast.cursor, 1)
            self.find_node(self.list_ast.cursor, model)
            return
        for ast in self.list_ast:
            output_cursor(ast.cursor, 1)
            self.find_node(ast.cursor, model)
        self.finalize(model)
        return model


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
        if namespace != self.parser_information["namespace"]:
            pass
        elif kind == CursorKind.ENUM_DECL and node.spelling != "": # alternative self.folder in node.location.file.name:
            model.enum_dict[node.spelling] = []
        elif kind == CursorKind.ENUM_CONSTANT_DECL and node.semantic_parent.spelling in model.enum_dict.keys():
            key = node.semantic_parent.spelling
            model.enum_dict[key].append(key + "::" + node.spelling)
        elif kind == CursorKind.CLASS_DECL:
            self.check_class(node, model)
        elif kind == CursorKind.CONSTRUCTOR:
            self.check_constructor(node, model)
        elif kind == CursorKind.STRUCT_DECL:
            self.check_struct(node, model)


    def check_class(self, node, model):
        """
        Inspect the nodes of kind CLASS_DECL and writes needed informations into model.
        Informations: model.name, model.population_groups
        """
        if node.spelling == self.parser_information["model_name"]:
            model.name = node.spelling
            base_classes = self.get_base_class_string(node)
            for base_class in base_classes:
                if base_class.startswith("CompartmentalModel:Populations:"):
                    model.population_groups.append(base_class.lstrip("CompartmentalModel:Populations:"))
        elif self.parser_information["optional"].get("simulation_name") and node.spelling == self.parser_information["optional"].get("simulation_name"):
            model.simulation_name = node.spelling
        
    
    def check_constructor(self, node, model):
        """
        Inspect the nodes of kind CONSTRUCTOR and writes needed informations into model.
        Informations: model.init
        """
        if node.spelling == model.name:
            init = {"type" : [], "name" : []}
            for arg in node.get_arguments():
                tokens = []
                for token in arg.get_tokens():
                    tokens.append(token.spelling)
                init["type"].append(" ".join(tokens[:-1]))
                init["name"].append(tokens[-1])
            model.init.append(init)

    def check_struct(self, node, model):
        if self.parser_information["parameterset"] in node.location.file.name:
            pass
    
    def finalize(self, model):
        for key, value in self.data_information.items():
            model.set_attribute(key, value)
        model.set_attribute("namespace", self.parser_information["namespace"])
        
        #assert(model.name != None), "set a model name"
        #assert(model.namespace != None), "set a model name"

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
        output_cursor_and_children(self.list_ast[index].cursor)
    
    def output_ast_file(self, index):
        with open('output.txt', 'a') as f:
            output_cursor_and_children_file(self.list_ast[index].cursor, f)

def indent(level):
    """ 
    Indentation string for pretty-printing
    """ 
    return '  '*level

def output_cursor(cursor, level):
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

    print(indent(level) + spelling, '<' + str(kind) + '>')
    print(indent(level+1) + '"'  + displayname + '"')


def output_cursor_and_children(cursor, level=0):
    """ 
    Output this cursor and its children with minimal formatting.
    """
    output_cursor(cursor, level)
    if cursor.kind.is_reference():
        print(indent(level) + 'reference to:')
        output_cursor(cursor.referenced, level+1)

    # Recurse for children of this cursor
    has_children = False
    for c in cursor.get_children():
        if not has_children:
            print(indent(level) + '{')
            has_children = True
        output_cursor_and_children(c, level+1)

    if has_children:
        print(indent(level) + '}')

def indent(level):
    """ 
    Indentation string for pretty-printing
    """ 
    return '  '*level

def output_cursor_file(cursor, f, level):
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

    f.write(indent(level) + spelling + ' <' + str(kind) + '> ')
    if cursor.location.file:
        f.write(cursor.location.file.name + '\n')
    f.write(indent(level+1) + '"'  + displayname + '"\n')

def output_cursor_and_children_file(cursor, f, level=0):
    """ 
    Output this cursor and its children with minimal formatting.
    """
    output_cursor_file(cursor, f, level)
    if cursor.kind.is_reference():
        f.write(indent(level) + 'reference to:\n')
        output_cursor_file(cursor.referenced, f, level+1)

    # Recurse for children of this cursor
    has_children = False
    for c in cursor.get_children():
        if not has_children:
            f.write(indent(level) + '{\n')
            has_children = True
        output_cursor_and_children_file(c, f, level+1)

    if has_children:
        f.write(indent(level) + '}\n')
