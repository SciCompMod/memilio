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