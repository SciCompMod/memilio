import os
from os.path import exists
from warnings import catch_warnings
from clang.cindex import *
import subprocess
from subprocess import check_output, call
from dataclasses_json import config
from memilio.generation import IntermediateRepresentation
import memilio.generation.utility as Utility
import tempfile

class Scanner:
    def __init__(self, conf):
        
        # Maybe use clang -v to get install dir from clang. 
        # Need to be same version as libclang
        self.config = conf
        Utility.try_set_libclang_path(self.config.optional.get("libclang_library_path"))
        self.ast = None
        self.create_ast()

    def create_ast(self):
        """
        Creates an AST for a main .cpp file with database. Requires an compile_commands.json.
        Saves them in list_ast.
        """
        idx = Index.create()
        
        file_args = []
        compdb = CompilationDatabase.fromDirectory(self.config.project_path + self.config.path_database)
        commands = compdb.getCompileCommands(self.config.project_path + self.config.source_file)
        for command in commands:
          for argument in command.arguments:
                if argument != '-Wno-unknown-warning' and argument!= '/usr/bin/g++' and argument != "--driver-mode=g++":
                    file_args.append(argument)
        file_args = file_args[:-4]

        clang_cmd = ["clang", self.config.project_path + self.config.source_file, '-emit-ast', '-o', '-']
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
        Extracts the information of the asts and saves them in the data class intermed_repr.
        Iterates over list of list_ast and calls find_node to visit all nodes of ast.
        """
        intermed_repr = IntermediateRepresentation()
        Utility.output_cursor(self.ast.cursor, 1)
        self.find_node(self.ast.cursor, intermed_repr)
        self.finalize(intermed_repr)
        return intermed_repr

    def find_node(self, node, intermed_repr, namespace = ""):
        """
        Recursively walks over every node of an ast. Saves the namespace the node is in.
        Calls check_node_kind for extracting information out of the nodes.
        """
        if node.kind == CursorKind.NAMESPACE:
            namespace = (namespace + node.spelling + "::")
        elif namespace == self.config.namespace:
            self.switch_node_kind(node.kind)(node, intermed_repr)
        
        for n in node.get_children(): 
            self.find_node(n, intermed_repr, namespace)

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
            CursorKind.STRUCT_DECL: self.check_struct,
            CursorKind.TYPE_ALIAS_DECL: self.check_type_alias
        }
        return switch.get(kind, lambda *args: None)

    def check_enum(self, node, intermed_repr):
        if node.spelling.strip() != "": # alternative self.folder in node.location.file.name:
            intermed_repr.enum_populations[node.spelling] = []

    def check_enum_const(self, node, intermed_repr):
        if node.semantic_parent.spelling in intermed_repr.enum_populations.keys():
            key = node.semantic_parent.spelling
            intermed_repr.enum_populations[key].append(node.spelling)

    def check_class(self, node, intermed_repr):
        """
        Inspect the nodes of kind CLASS_DECL and writes needed informations into intermed_repr.
        Informations: intermed_repr.model_class, intermed_repr.population_groups
        """
        if node.spelling == self.config.model_class:
            intermed_repr.model_class = node.spelling
            self.check_model_base(node, intermed_repr)
            if self.config.optional.get("age_group"):
                self.check_age_group(node, intermed_repr)
        elif self.config.optional.get("simulation_class") and node.spelling == self.config.optional.get("simulation_class"):
            intermed_repr.simulation_class = node.spelling
        elif self.config.optional.get("parameterset_wrapper") and self.config.namespace + self.config.parameterset in [base.spelling for base in node.get_children()]:
            intermed_repr.parameterset_wrapper = node.spelling
    
    def check_model_base(self, node, intermed_repr):

        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            base_type = base.get_definition().type
            intermed_repr.model_base = Utility.get_base_class_string(base_type)
    
    def check_base_specifier(self, node, intermed_repr):
        pass

    def check_age_group(self, node, intermed_repr):
        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            for base_template_arg in base.get_children():
                if base_template_arg.kind == CursorKind.TYPE_REF and "AgeGroup" in base_template_arg.spelling:
                    for child in base_template_arg.get_definition().get_children():
                        if child.kind == CursorKind.CXX_BASE_SPECIFIER:
                            intermed_repr.age_group["base"] = child.get_definition().type.spelling
                        elif child.kind == CursorKind.CONSTRUCTOR:
                            intermed_repr.age_group["init"] = [arg.spelling for arg in child.type.argument_types()]
                            

    def check_constructor(self, node, intermed_repr):
        """
        Inspect the nodes of kind CONSTRUCTOR and writes needed informations into intermed_repr.
        Informations: intermed_repr.init
        """
        if node.spelling == intermed_repr.model_class:
            init = {"type" : [], "name" : []}
            for arg in node.get_arguments():
                tokens = []
                for token in arg.get_tokens():
                    tokens.append(token.spelling)
                init["type"].append(" ".join(tokens[:-1]))
                init["name"].append(tokens[-1])
            intermed_repr.model_init.append(init)

    def check_type_alias(self, node, intermed_repr):
        if node.spelling == self.config.parameterset:
            intermed_repr.parameterset = node.spelling
    
    def check_struct(self, node, intermed_repr):
        pass
    
    def finalize(self, intermed_repr):

        # remove unnecesary enum
        population_groups = []
        for value in intermed_repr.model_base[1:]:
            if "Population" in value[0]:
                population_groups = [pop[0].split("::")[-1] for pop in value[1:]]
        intermed_repr.population_groups = population_groups
        new_enum = {}
        for key in intermed_repr.enum_populations:
            if key in population_groups:
                new_enum[key] = intermed_repr.enum_populations[key]
        intermed_repr.enum_populations = new_enum

        intermed_repr.set_attribute("namespace", self.config.namespace)
        intermed_repr.set_attribute("python_module_name", self.config.optional.get("python_module_name"))
        intermed_repr.set_attribute("target_folder", self.config.target_folder)
        intermed_repr.set_attribute("project_path", self.config.project_path)
        
        intermed_repr.check_complete_data(self.config.optional)

    def output_ast(self):
        """
        Outputs the ast
        """
        Utility.output_cursor_and_children(self.ast.cursor)
    
    def output_ast_file(self):
        with open('output_ast.txt', 'a') as f:
            Utility.output_cursor_and_children_file(self.ast.cursor, f)
