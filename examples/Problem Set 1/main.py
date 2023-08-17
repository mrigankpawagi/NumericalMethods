from solution import Solution

answers = Solution.solve_all(10)
for i, ans in enumerate(answers):
    print(f'Problem {i + 1}\t {ans}')