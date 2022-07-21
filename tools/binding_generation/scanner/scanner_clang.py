import os
from os.path import exists
from clang.cindex import *
import subprocess
from subprocess import check_output, call
import clang.cindex
from dataclasses_json import config
from generator.Model import Model
from scanner.utility import *
import tempfile
import re

class Scanner:
    def __init__(self, conf):
        
        # Maybe use clang -v to get install dir from clang. Need to be same version
        self.config = conf
        Config.set_library_file(self.config.libclang_library_path)

        self.ast = None
        self.create_ast()
        for dig in self.ast.diagnostics:
            print(dig)
        #output_cursor_and_children(list_ast[1].cursor)

    def create_ast(self):
        """
        Creates an Ast for a main .cpp file with database. Requires an compile_commands.json.
        Saves them in list_ast.
        """
        # look if folder exists
        idx = Index.create()
        
        file_args = []
        compdb = CompilationDatabase.fromDirectory(check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + self.config.path_database)
        commands = compdb.getCompileCommands(self.config.source_file)
        #file_args = ["clang++"]#, folder_path + "/" + file_name]
        for command in commands:
          for argument in command.arguments:
                if argument != '-Wno-unknown-warning' and argument!= '/usr/bin/g++' and argument != "--driver-mode=g++":
                    file_args.append(argument)
        file_args = file_args[:-4]

        clang_cmd = ["clang", self.config.source_file, '-emit-ast', '-o', '-']
        clang_cmd.extend(file_args)
        clang_cmd_result = subprocess.run(clang_cmd, stdout=subprocess.PIPE)
        clang_cmd_result.check_returncode()
        with tempfile.NamedTemporaryFile() as ast_file:
            # Since `clang.Index.read` expects a file path, write generated AST to a
            # temporary named file. This file will be automatically deleted when closed.
            ast_file.write(clang_cmd_result.stdout)
            self.ast = idx.read(ast_file.name)

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
        elif namespace == self.config.namespace:
            self.switch_node_kind(node.kind)(node, model)
        
        for n in node.get_children(): 
            self.find_node(n, model, namespace)

    def switch_node_kind(self, kind):
        """
        Returns the appropriate method for the given kind.
        """
        switch = {
            CursorKind.ENUM_DECL: self.check_enum,
            CursorKind.ENUM_CONSTANT_DECL: self.check_enum_const,
            CursorKind.CLASS_DECL: self.check_class,
            CursorKind.CLASS_TEMPLATE: self.check_class,
            CursorKind.CXX_BASE_SPECIFIER: self.check_base_specifier,
            CursorKind.CONSTRUCTOR: self.check_constructor,
            CursorKind.STRUCT_DECL: self.check_struct
        }
        return switch.get(kind, lambda *args: None)

    def check_enum(self, node, model):
        if node.spelling.strip() != "": # alternative self.folder in node.location.file.name:
            model.enum_populations[node.spelling] = []

    def check_enum_const(self, node, model):
        if node.semantic_parent.spelling in model.enum_populations.keys():
            key = node.semantic_parent.spelling
            model.enum_populations[key].append(key + "::" + node.spelling)

    def check_class(self, node, model):
        """
        Inspect the nodes of kind CLASS_DECL and writes needed informations into model.
        Informations: model.name, model.population_groups
        """
        if node.spelling == self.config.model_class:
            model.model_class = node.spelling
            self.check_base_specifier(node, model)
            if self.config.optional.get("age_group"):
                self.check_age_group(node, model)
        elif self.config.optional.get("simulation_class") and node.spelling == self.config.optional.get("simulation_class"):
            model.simulation_class = node.spelling
    
    def check_base_specifier(self, node, model):

        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            base_type = base.get_definition().type
            model.model_base = get_base_class_string(base_type)
    
    def check_age_group(self, node, model):

        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            for base_template_arg in base.get_children():
                if base_template_arg.kind == CursorKind.TYPE_REF and "AgeGroup" in base_template_arg.spelling:
                    for child in base_template_arg.get_definition().get_children():
                        if child.kind == CursorKind.CXX_BASE_SPECIFIER:
                            model.age_group["base"] = child.get_definition().type.spelling
                        elif child.kind == CursorKind.CONSTRUCTOR:
                            model.age_group["init"] = [arg.spelling for arg in child.type.argument_types()]
                            

    def check_constructor(self, node, model):
        """
        Inspect the nodes of kind CONSTRUCTOR and writes needed informations into model.
        Informations: model.init
        """
        if node.spelling == model.model_class:
            init = {"type" : [], "name" : []}
            for arg in node.get_arguments():
                tokens = []
                for token in arg.get_tokens():
                    tokens.append(token.spelling)
                init["type"].append(" ".join(tokens[:-1]))
                init["name"].append(tokens[-1])
            model.model_init.append(init)

    def check_struct(self, node, model):
        #if self.config.parameterset_name in node.location.file.name:
        pass
    
    def finalize(self, model):

        # remove unnecesary enum
        population_groups = []
        for value in model.model_base[1:]:
            if "Population" in value[0]:
                population_groups = [pop[0].split("::")[-1] for pop in value[1:]]
        model.population_groups = population_groups
        new_enum = {}
        for key in model.enum_populations:
            if key in population_groups:
                new_enum[key] = model.enum_populations[key]
        model.enum_populations = new_enum

        model.set_attribute("namespace", self.config.namespace)
        model.set_attribute("python_module_name", self.config.optional.get("python_module_name"))
        
        #assert(model.name != None), "set a model name"
        #assert(model.namespace != None), "set a model name"

    def output_ast(self):
        """
        Outputs the ast
        """
        output_cursor_and_children(self.ast.cursor)
    
    def output_ast_file(self):
        with open('output_ast.txt', 'a') as f:
            output_cursor_and_children_file(self.ast.cursor, f)