from __future__ import annotations

from collections import defaultdict, deque
from typing import Iterable

from .errors import ValidationError
from .models import PedigreeIndividual


class Pedigree:
    """Pedigree relationship and inbreeding calculations.

    The numerator relationship matrix A is built in topological order. The
    expected inbreeding coefficient of a future offspring is A[sire, dam] / 2.
    Unknown founders are treated as unrelated in this prototype.
    """

    def __init__(self, individuals: Iterable[PedigreeIndividual]):
        entries = list(individuals)
        if not entries:
            raise ValidationError("pedigree_individuals_required")
        self._individuals = {entry.individual_id: entry for entry in entries}
        if len(self._individuals) != len(entries):
            raise ValidationError("pedigree_duplicate_individual_id")
        self._order = self._topological_order()
        self._index = {individual_id: index for index, individual_id in enumerate(self._order)}
        self._matrix = self._build_relationship_matrix()

    @property
    def order(self) -> list[str]:
        return list(self._order)

    def _topological_order(self) -> list[str]:
        indegree: dict[str, int] = {individual_id: 0 for individual_id in self._individuals}
        children: dict[str, list[str]] = defaultdict(list)

        for entry in self._individuals.values():
            for parent_id in (entry.sire_id, entry.dam_id):
                if parent_id is None or parent_id not in self._individuals:
                    continue
                indegree[entry.individual_id] += 1
                children[parent_id].append(entry.individual_id)

        queue = deque(
            individual_id for individual_id, degree in indegree.items() if degree == 0
        )
        order: list[str] = []
        while queue:
            individual_id = queue.popleft()
            order.append(individual_id)
            for child_id in children[individual_id]:
                indegree[child_id] -= 1
                if indegree[child_id] == 0:
                    queue.append(child_id)

        if len(order) != len(self._individuals):
            raise ValidationError("pedigree_cycle_detected")
        return order

    def _build_relationship_matrix(self) -> list[list[float]]:
        size = len(self._order)
        matrix = [[0.0 for _ in range(size)] for _ in range(size)]

        for i, individual_id in enumerate(self._order):
            entry = self._individuals[individual_id]
            sire_index = self._index.get(entry.sire_id) if entry.sire_id else None
            dam_index = self._index.get(entry.dam_id) if entry.dam_id else None

            for j in range(i):
                sire_value = matrix[sire_index][j] if sire_index is not None else 0.0
                dam_value = matrix[dam_index][j] if dam_index is not None else 0.0
                matrix[i][j] = matrix[j][i] = 0.5 * (sire_value + dam_value)

            if sire_index is not None and dam_index is not None:
                matrix[i][i] = 1.0 + 0.5 * matrix[sire_index][dam_index]
            else:
                matrix[i][i] = 1.0

        return matrix

    def relationship(self, first_id: str, second_id: str) -> float:
        try:
            return self._matrix[self._index[first_id]][self._index[second_id]]
        except KeyError as exc:
            raise ValidationError(f"pedigree_individual_unknown:{exc.args[0]}") from exc

    def kinship(self, first_id: str, second_id: str) -> float:
        return 0.5 * self.relationship(first_id, second_id)

    def individual_inbreeding(self, individual_id: str) -> float:
        index = self._index.get(individual_id)
        if index is None:
            raise ValidationError(f"pedigree_individual_unknown:{individual_id}")
        return self._matrix[index][index] - 1.0

    def offspring_inbreeding(self, sire_id: str, dam_id: str) -> float:
        return 0.5 * self.relationship(sire_id, dam_id)

    def ancestors(self, individual_id: str) -> set[str]:
        if individual_id not in self._individuals:
            raise ValidationError(f"pedigree_individual_unknown:{individual_id}")
        result: set[str] = set()
        stack = [individual_id]
        while stack:
            current = stack.pop()
            entry = self._individuals[current]
            for parent_id in (entry.sire_id, entry.dam_id):
                if parent_id and parent_id in self._individuals and parent_id not in result:
                    result.add(parent_id)
                    stack.append(parent_id)
        return result

    def is_ancestor(self, candidate_ancestor_id: str, descendant_id: str) -> bool:
        return candidate_ancestor_id in self.ancestors(descendant_id)

    def are_full_siblings(self, first_id: str, second_id: str) -> bool:
        if first_id == second_id:
            return False
        try:
            first = self._individuals[first_id]
            second = self._individuals[second_id]
        except KeyError as exc:
            raise ValidationError(f"pedigree_individual_unknown:{exc.args[0]}") from exc
        return (
            first.sire_id is not None
            and first.dam_id is not None
            and first.sire_id == second.sire_id
            and first.dam_id == second.dam_id
        )

    def hard_relation_block(self, male_id: str, female_id: str) -> str | None:
        if self.is_ancestor(male_id, female_id) or self.is_ancestor(female_id, male_id):
            return "direct_ancestor_pair"
        if self.are_full_siblings(male_id, female_id):
            return "full_sibling_pair"
        return None

    def matrix_as_dict(self) -> dict[str, dict[str, float]]:
        return {
            row_id: {
                column_id: self._matrix[row_index][column_index]
                for column_index, column_id in enumerate(self._order)
            }
            for row_index, row_id in enumerate(self._order)
        }
