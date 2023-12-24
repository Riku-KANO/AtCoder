from dataclasses import dataclass
from typing import List, Tuple
import sys
import numpy as np

MAX_INVEST_LEVEL = 20

INVEST = 4

@dataclass
class Project:
    h: int
    v: int


@dataclass
class Card:
    t: int
    w: int
    p: int


class Judge:

    def __init__(self, n: int, m: int, k: int):
        self.n = n
        self.m = m
        self.k = k

    def read_initial_cards(self) -> List[Card]:
        cards = []
        for _ in range(self.n):
            t, w = map(int, input().split())
            cards.append(Card(t, w, 0))
        return cards

    def read_projects(self) -> List[Project]:
        projects = []
        for _ in range(self.m):
            h, v = map(int, input().split())
            projects.append(Project(h, v))
        return projects

    def use_card(self, c: int, m: int) -> None:
        print(f"{c} {m}", flush=True)

    def read_money(self) -> int:
        return int(input())

    def read_next_cards(self) -> List[Card]:
        cards = []
        for _ in range(self.k):
            t, w, p = map(int, input().split())
            cards.append(Card(t, w, p))
        return cards

    def select_card(self, r: int) -> None:
        print(r, flush=True)

    def comment(self, message: str) -> None:
        print(f"# {message}")


class Solver:

    def __init__(self, n: int, m: int, k: int, t: int):
        self.n = n
        self.m = m
        self.k = k
        self.t = t
        self.judge = Judge(n, m, k)
        self.types = [0, 1, 2, 3, 4]
        self.hiddens = self._generate_hiddens()

    def solve(self) -> int:
        self.turn = 0
        self.money = 0
        self.invest_level = 0
        self.cards = self.judge.read_initial_cards()
        self.projects = self.judge.read_projects()


        for _ in range(self.t):
            use_card_i, use_target = self._select_action()
            if self.cards[use_card_i].t == INVEST:
                self.invest_level += 1
            # example for comments
            self.judge.comment(f"used {self.cards[use_card_i]} to target {use_target}")
            self.judge.use_card(use_card_i, use_target)
            assert self.invest_level <= MAX_INVEST_LEVEL

            self.projects = self.judge.read_projects()
            self.money = self.judge.read_money()

            next_cards = self.judge.read_next_cards()
            select_card_i = self._select_next_card(next_cards)
            self.cards[use_card_i] = next_cards[select_card_i]
            self.judge.select_card(select_card_i)
            self.money -= next_cards[select_card_i].p
            assert self.money >= 0

            self.turn += 1

        return self.money

    def _select_action(self) -> Tuple[int, int]:
        # TODO: implement your strategy
        project_values = [project.v / project.h / 1.25 for project in self.projects]
        values = []
        for idx, card in enumerate(self.cards):
            if card.t == 0:
                values.append(max(project_values))
            elif card.t == 2:
                values.append(1 / min(project_values))
            elif card.t == 1:
                values.append(sum(map(lambda x: x, project_values))/len(project_values))
            elif card.t == 3:
                values.append(sum(map(lambda x: 1/x, project_values))/len(project_values))
            else:
                values.append(1.5)
        card_value_pairs = [(idx,val) for idx,val in enumerate(values)]
        card_value_pairs.sort(key = lambda x: x[1])
        target_card_idx = card_value_pairs[-1][0]
        target_card = self.cards[target_card_idx]
        target_project = 0  
        if target_card.t == 0 or target_card.t == 2:
            target_project = np.argmax(project_values)

        return (target_card_idx, target_project)

    def _select_next_card(self, next_cards: List[Card]) -> int:
        # TODO: implement your strategy
        project_values = [project.v / project.h / 1.25 for project in self.projects]
        values = []
        for idx, card in enumerate(next_cards):
            
            if card.p > self.money:
                values.append(-100)
                continue

            if card.t == 0 or card.t == 1:
                denom = card.p if card.p != 0 else 1
                values.append(card.w / denom * 0.8)
            elif card.t == 2:
                values.append(1 / min(project_values))
            elif card.t == 3:
                values.append(sum(map(lambda x: 1/x, project_values))/len(project_values))
            elif card.t == 4:
                values.append(1.5)
        
        return np.argmax(values)
    
    def _generate_hiddens(self, size=1):
        x0 = np.random.randint(1,20, size=size)
        x1 = np.random.randint(1,10, size=size)
        x2 = np.random.randint(1,10, size=size)
        x3 = np.random.randint(1,5, size=size)
        x4 = np.random.randint(1,3, size=size)
        if size == 1:
            return np.array([x0, x1, x2, x3, x4])
        else:
            return np.array([list(tup) for tup in zip(x0,x1,x2,x3,x4)])
        
    def _get_available_card_idx(self, cards: List[Card]) -> List[int]:
        return [idx for idx, card in enumerate(cards) if card.p <= self.money]



def main():
    n, m, k, t = map(int, input().split())
    solver = Solver(n, m, k, t)
    score = solver.solve()
    print(f"score:{score}", file=sys.stderr)


if __name__ == "__main__":
    main()
