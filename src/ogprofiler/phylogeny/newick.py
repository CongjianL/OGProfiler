"""Small dependency-free Newick tree model, parser, pruning, and rerooting."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from ogprofiler.exceptions import PhylogenyError


@dataclass(slots=True)
class TreeNode:
    name: str | None = None
    length: float | None = None
    children: list[TreeNode] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


class _Parser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.position = 0

    def parse(self) -> TreeNode:
        if not self.text:
            raise PhylogenyError("Newick input is empty")
        root = self.node()
        self.space()
        if self.position >= len(self.text) or self.text[self.position] != ";":
            raise PhylogenyError("Newick tree must end with ';'")
        self.position += 1
        self.space()
        if self.position != len(self.text):
            raise PhylogenyError("Unexpected content after Newick tree")
        return root

    def space(self) -> None:
        while self.position < len(self.text) and self.text[self.position].isspace():
            self.position += 1

    def node(self) -> TreeNode:
        self.space()
        children: list[TreeNode] = []
        if self.position < len(self.text) and self.text[self.position] == "(":
            self.position += 1
            while True:
                children.append(self.node())
                self.space()
                if self.position >= len(self.text):
                    raise PhylogenyError("Unclosed Newick child list")
                token = self.text[self.position]
                self.position += 1
                if token == ")":
                    break
                if token != ",":
                    raise PhylogenyError(f"Unexpected Newick token: {token}")
        name = self.label()
        length = None
        self.space()
        if self.position < len(self.text) and self.text[self.position] == ":":
            self.position += 1
            raw = self.unquoted(stop=",();")
            try:
                length = float(raw)
            except ValueError as error:
                raise PhylogenyError(f"Invalid Newick branch length: {raw}") from error
        if not children and not name:
            raise PhylogenyError("Newick leaf is missing a name")
        return TreeNode(name or None, length, children)

    def label(self) -> str:
        self.space()
        if self.position >= len(self.text):
            return ""
        if self.text[self.position] == "'":
            self.position += 1
            value: list[str] = []
            while self.position < len(self.text):
                token = self.text[self.position]
                self.position += 1
                if token != "'":
                    value.append(token)
                    continue
                if self.position < len(self.text) and self.text[self.position] == "'":
                    value.append("'")
                    self.position += 1
                    continue
                return "".join(value)
            raise PhylogenyError("Unclosed quoted Newick label")
        return self.unquoted(stop=":,();")

    def unquoted(self, stop: str) -> str:
        self.space()
        start = self.position
        while self.position < len(self.text) and self.text[self.position] not in stop:
            self.position += 1
        return self.text[start : self.position].strip()


def parse_newick(text: str) -> TreeNode:
    return _Parser(text).parse()


def leaf_names(root: TreeNode) -> tuple[str, ...]:
    names: list[str] = []

    def visit(node: TreeNode) -> None:
        if node.is_leaf:
            if node.name is None:
                raise PhylogenyError("Tree leaf is missing a name")
            names.append(node.name)
            return
        for child in node.children:
            visit(child)

    visit(root)
    if len(names) != len(set(names)):
        raise PhylogenyError("Tree contains duplicate leaf names")
    return tuple(names)


def _quoted(name: str) -> str:
    if any(token in name for token in "(),:;[]' \t\r\n"):
        return "'" + name.replace("'", "''") + "'"
    return name


def to_newick(root: TreeNode) -> str:
    def render(node: TreeNode) -> str:
        prefix = (
            "(" + ",".join(render(child) for child in node.children) + ")"
            if node.children
            else ""
        )
        label = _quoted(node.name) if node.name else ""
        length = f":{node.length:.12g}" if node.length is not None else ""
        return prefix + label + length

    return render(root) + ";"


def prune_tree(root: TreeNode, keep: set[str]) -> TreeNode:
    def prune(node: TreeNode) -> TreeNode | None:
        if node.is_leaf:
            return TreeNode(node.name, node.length) if node.name in keep else None
        children = [value for child in node.children if (value := prune(child)) is not None]
        if not children:
            return None
        if len(children) == 1:
            child = children[0]
            if node.length is not None:
                child.length = (child.length or 0.0) + node.length
            return child
        return TreeNode(node.name, node.length, children)

    result = prune(root)
    if result is None:
        raise PhylogenyError("Species-tree pruning removed every leaf")
    result.length = None
    return result


def _graph(root: TreeNode) -> tuple[list[TreeNode], dict[int, list[tuple[int, float]]]]:
    nodes: list[TreeNode] = []
    adjacency: dict[int, list[tuple[int, float]]] = {}

    def visit(node: TreeNode, parent: int | None = None) -> int:
        node_id = len(nodes)
        nodes.append(node)
        adjacency[node_id] = []
        if parent is not None:
            length = node.length if node.length is not None else 1.0
            adjacency[node_id].append((parent, length))
            adjacency[parent].append((node_id, length))
        for child in node.children:
            visit(child, node_id)
        return node_id

    visit(root)
    return nodes, adjacency


def _orient(
    nodes: list[TreeNode],
    adjacency: dict[int, list[tuple[int, float]]],
    root_edge: tuple[int, int, float],
) -> TreeNode:
    left, right, left_length = root_edge
    edge_length = next(length for neighbor, length in adjacency[left] if neighbor == right)
    right_length = edge_length - left_length

    def clone(node_id: int, parent: int, length: float) -> TreeNode:
        source = nodes[node_id]
        children = [
            clone(neighbor, node_id, branch)
            for neighbor, branch in adjacency[node_id]
            if neighbor != parent
        ]
        return TreeNode(source.name, length, children)

    return TreeNode(
        children=[clone(left, right, left_length), clone(right, left, right_length)]
    )


def midpoint_root(root: TreeNode) -> TreeNode:
    nodes, adjacency = _graph(root)
    leaves = [node_id for node_id, neighbors in adjacency.items() if len(neighbors) == 1]
    if len(leaves) < 2:
        return root

    def distances(start: int) -> tuple[dict[int, float], dict[int, int]]:
        values = {start: 0.0}
        parents: dict[int, int] = {}
        stack = [start]
        while stack:
            node_id = stack.pop()
            for neighbor, length in adjacency[node_id]:
                if neighbor in values:
                    continue
                values[neighbor] = values[node_id] + length
                parents[neighbor] = node_id
                stack.append(neighbor)
        return values, parents

    first_distances, _ = distances(leaves[0])
    first = max(leaves, key=lambda value: first_distances[value])
    values, parents = distances(first)
    second = max(leaves, key=lambda value: values[value])
    path = [second]
    while path[-1] != first:
        path.append(parents[path[-1]])
    target = values[second] / 2.0
    travelled = 0.0
    for left, right in zip(path, path[1:], strict=False):
        length = next(value for neighbor, value in adjacency[left] if neighbor == right)
        if travelled + length >= target:
            offset = target - travelled
            return _orient(nodes, adjacency, (left, right, offset))
        travelled += length
    raise PhylogenyError("Failed to locate gene-tree midpoint")


def outgroup_root(root: TreeNode, outgroup: str) -> TreeNode:
    nodes, adjacency = _graph(root)
    matches = [index for index, node in enumerate(nodes) if node.is_leaf and node.name == outgroup]
    if len(matches) != 1:
        raise PhylogenyError(f"Outgroup leaf not found uniquely: {outgroup}")
    leaf = matches[0]
    neighbor, length = adjacency[leaf][0]
    return _orient(nodes, adjacency, (leaf, neighbor, length / 2.0))


def edge_rootings(root: TreeNode) -> Iterator[TreeNode]:
    nodes, adjacency = _graph(root)
    for left in sorted(adjacency):
        for right, length in adjacency[left]:
            if left < right:
                yield _orient(nodes, adjacency, (left, right, length / 2.0))
