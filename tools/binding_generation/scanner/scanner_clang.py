import os
from os.path import exists
from clang.cindex import *
import subprocess
from subprocess import check_output, call
import clang.cindex
from generator.Model import Model
from scanner.utility import *
import tempfile
import re

class Scanner:
    def __init__(self, conf):

        Config.set_library_file(conf["scanner_input"]["libclang_library_path"])
        self.folder = conf["scanner_input"]["source_folder"]
        self.parser_information = conf["parser_information"]
        self.data_information = conf["data_information"]
        self.database = conf["database"]
        
        self.ast = None
        self.create_ast(self.database["main_file"])
        for dig in self.ast.diagnostics:
            print(dig)
        #output_cursor_and_children(list_ast[1].cursor)

    def create_ast(self, file_name):
        """
        Creates an Ast for a main .cpp file with database. Requires an compile_commands.json.
        Saves them in list_ast.
        """
        # look if folder exists
        project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/models/"
        folder_path = project_path + self.folder
        idx = Index.create()
        
        file_args = []
        compdb = CompilationDatabase.fromDirectory(check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + self.database["path_database"])
        commands = compdb.getCompileCommands(folder_path + "/" + file_name)
        #file_args = ["clang++"]#, folder_path + "/" + file_name]
        for command in commands:
          for argument in command.arguments:
                if argument != '-Wno-unknown-warning' and argument!= '/usr/bin/g++' and argument != "--driver-mode=g++":
                    file_args.append(argument)
        file_args = file_args[:-4]

        clang_cmd = ["clang", folder_path + "/" + file_name, '-emit-ast', '-o', '-']
        clang_cmd.extend(file_args)
        clang_cmd_result = subprocess.run(clang_cmd, stdout=subprocess.PIPE)
        clang_cmd_result.check_returncode()
        with open('text.txt', 'wb', 0) as ast_file:
            # Since `clang.Index.read` expects a file path, write generated AST to a
            # temporary named file. This file will be automatically deleted when closed.
            ast_file.write(clang_cmd_result.stdout)
            self.ast = idx.read(ast_file.name)
        """
        clang_cmd_result = subprocess.Popen(["clang", "-v", "-x", "c++", "-"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, out = clang_cmd_result.communicate('')
        reg = re.compile('lib/clang.*/include$')
        print(out)
        #c_command = next(line.strip() for line in out.split('\n') if reg.search(line))
        print("finished")
        # Step 1: load the compilation database
        compdb = CompilationDatabase.fromDirectory(check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + self.database["path_database"])
        commands = compdb.getCompileCommands(folder_path + "/" + file_name)
        file_args = ["clang", "-cc1"]
        #file_args = ["clang++"]#, folder_path + "/" + file_name]
        for command in commands:
          for argument in command.arguments:
                if argument != '-Wno-unknown-warning' and argument!= '/usr/bin/g++':# and argument != "--driver-mode=g++":
                    file_args.append(argument)

        file_args = file_args[:-1]
        #file_args.append(r"-stdlib=libc++")
        #file_args.append(r'-nostdinc++')
        #file_args.append(r'/home/betz_mx/memilio/tools/binding_generation/out.txt')
        #file_args.append(r"-emit-ast")
        #file_args = file_args + c_commands
        file_args.append(r'-Xclang')
        file_args.append(r"-cc1")
        file_args.append(r'-frecovery-ast')
        #file_args.append(r"-stdlib=libc++")
        #file_args.append(r'-ast-dump')
        file_args.append(r'-Wdeprecated')
        file_args.append(r'-fsyntax-only')
        #file_args.append(r'-Wunused-command-line-argument')
        #call(file_args)
        #with open('output_command.txt', 'wb', 0) as outputfile:
        #    call(file_args, stdout=outputfile)
        self.ast = idx.parse(folder_path + "/" + file_name, file_args)
        self.ast.save("text.c")
        self.ast = idx.read("text.c")
        """
        
    def extract_results(self):
        """
        Extracts the information of the asts and saves them in the data class model.
        Iterates over list of list_ast and calls find_node to visit all nodes of ast.
        """
        model = Model()
        output_cursor(self.ast.cursor, 1)
        self.find_node(self.ast.cursor, model)
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

    def output_ast(self):
        """
        Outputs the ast
        """
        output_cursor_and_children(self.ast.cursor)
    
    def output_ast_file(self):
        with open('output_ast.txt', 'a') as f:
            output_cursor_and_children_file(self.ast.cursor, f)
