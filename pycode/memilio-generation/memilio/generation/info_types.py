from dataclasses import dataclass, field


@dataclass
class method_type_info:
    type: str
    name: str
    cursorkind: str
    namespace: str
    return_type: str
    arg_types: list[str] = field(default_factory=list)
    arg_names: list[str] = field(default_factory=list)
    parent_name: str = ""
    is_const: bool = False
    is_member: bool = False

    def signature_key(self) -> tuple:
        return (
            self.name,
            self.cursorkind,
            tuple(self.arg_types),
            tuple(self.arg_names),
            self.parent_name,
        )


@dataclass
class binding_type_info:
    type: str
    name: str
    cursorkind: str
    namespace: str
    return_type: str
    arg_types: list[str] = field(default_factory=list)
    arg_names: list[str] = field(default_factory=list)
    parent_name: str = ""
    is_const: bool = False
    is_member: bool = False
    methods: list[method_type_info] = field(default_factory=list)
    base_classes: list[str] = field(default_factory=list)
    init: list[dict[str, str]] = field(default_factory=list)
    template_args: list[str] = field(default_factory=list)

    def signature_key(self) -> tuple:
        return (
            self.name,
            self.cursorkind,
            tuple(self.arg_types),
            tuple(self.arg_names),
            self.parent_name,
        )
