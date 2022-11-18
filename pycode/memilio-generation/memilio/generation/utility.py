#############################################################################
# Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
#
# Authors: Maximilian Betz
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""
@file utility.py
@brief Additional functions used for the code generation.
"""
from typing import Any, List, TextIO
from clang.cindex import Config, Cursor, Type
import subprocess
import os


def try_set_libclang_path(path: str) -> None:
    """
    Try to set the file_path for the libclang library. 
    If its already set, the returned Exception gets caught and discarded.
    If the given path string is empty, it is obtain with a call on the terminal.

    @param path Path to the library files of libClang. Can be an empty string.
    """
    # Check if path was set in config. If not, try to get it with cmd.
    if (not path):
        clang_cmd = ["clang", '-print-file-name=libclang.so']
        clang_cmd_result = subprocess.check_output(clang_cmd)
        path = clang_cmd_result.rstrip()
    try:
        Config.set_library_file(os.path.abspath(path))
    except Exception as e:
        if str(e) != "library file must be set before before using any other functionalities in libclang.":
            raise (e)


def get_base_class(base_type: Type) -> List[Any]:
    """
    Retrieve the base class.
    Example for base_type: CompartmentalModel.

    @param Type of the current node.
    """
    result = [base_type.replace('> >', '>>')]
    for i in range(base_type.get_num_template_arguments()):
        result.append(get_base_class(base_type.get_template_argument_type(i)))
    return result


def get_base_class_string(base_type: Type) -> List[Any]:
    """
    Retrieve the spelling of the base class.
    Example for base_type.spelling: CompartmentalModel<mio::Populations<mio::AgeGroup, mio::InfectionState>, mio::SecirParams>.

    @param Type of the current node.
    """
    result = [base_type.spelling.replace('> >', '>>')]
    for i in range(base_type.get_num_template_arguments()):
        result.append(get_base_class_string(
            base_type.get_template_argument_type(i)))
    return result


def indent(level: int) -> str:
    """ 
    Indentation string for pretty-printing.

    @param level Amount of indentations.
    """
    return '  '*level


def output_cursor_print(cursor: Cursor, level: int) -> None:
    """ 
    Low level cursor output to the terminal.

    @param cursor Represents the current node of the AST as an Cursor object from libClang.
    @param level Amount of indentations.
    """
    spelling = ''
    displayname = ''

    if cursor.spelling:
        spelling = cursor.spelling
    if cursor.displayname:
        displayname = cursor.displayname
    kind = cursor.kind

    print(indent(level) + spelling, '<' + str(kind) + '>')
    print(indent(level+1) + '"' + displayname + '"')


def output_cursor_and_children_print(cursor: Cursor, level: int = 0) -> None:
    """ 
    Output this cursor and its children with minimal formatting to the terminal.

    @param cursor Represents the current node of the AST as an Cursor object from libClang.
    @param level [Default = 0] Amount of indentations.
    """
    output_cursor_print(cursor, level)
    if cursor.kind.is_reference():
        print(indent(level) + 'reference to:')
        output_cursor_print(cursor.referenced, level+1)

    # Recurse for children of this cursor
    has_children = False
    for c in cursor.get_children():
        if not has_children:
            print(indent(level) + '{')
            has_children = True
        output_cursor_and_children_print(c, level+1)

    if has_children:
        print(indent(level) + '}')


def output_cursor_file(cursor: Cursor, f: TextIO, level: int) -> None:
    """ 
    Low level cursor output to a file.

    @param cursor Represents the current node of the AST as an Cursor object from libClang.
    @param f Open file object for output.
    @param level Amount of indentations.
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
    f.write(indent(level+1) + '"' + displayname + '"\n')


def output_cursor_and_children_file(
        cursor: Cursor, f: TextIO, level: int = 0) -> None:
    """ 
    Output this cursor and its children with minimal formatting to a file.

    @param cursor Represents the current node of the AST as an Cursor object from libClang.
    @param f Open file object for output.
    @param level Amount of indentations.
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
