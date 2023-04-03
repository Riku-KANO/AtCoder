class RollingHash2D:
    def __init__(self, matrix):
        self.matrix = matrix
        self.base1 = 1007
        self.base2 = 2009
        self.mod1 = 10 ** 9 + 7
        self.mod2 = 10 ** 9 + 9
        self.h1 = [[0] * (len(matrix[0]) + 1) for _ in range(len(matrix))]
        self.h2 = [[0] * (len(matrix[0]) + 1) for _ in range(len(matrix))]
        self.p1 = [1] * (len(matrix[0]) + 1)
        self.p2 = [1] * (len(matrix) + 1)

        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                self.h1[i].append((self.h1[i][-1] * self.base1 + ord(matrix[i][j])) % self.mod1)
        for j in range(len(matrix[0])):
            for i in range(len(matrix)):
                self.h2[i].append((self.h2[i - 1][j] * self.base2 + self.h1[i][j]) % self.mod2)
        for i in range(1, len(matrix[0]) + 1):
            self.p1[i] = self.p1[i - 1] * self.base1 % self.mod1
        for i in range(1, len(matrix) + 1):
            self.p2[i] = self.p2[i - 1] * self.base2 % self.mod2

    def get_hash(self, x1, y1, x2, y2):
        h = (self.h2[x2][y2] - self.h2[x1 - 1][y2] - self.h2[x2][y1 - 1] + self.h2[x1 - 1][y1 - 1]) % self.mod2
        return h

if __name__ == "__main__":
  matrix1 = [
    "abcde",
    "fghij",
    "klmno",
    "pqrst",
    "uvwxy"
  ]

  matrix2 = [
      "abcde",
      "fghij",
      "klmno",
      "pqrst",
      "uvwxy"
  ]

  rh1 = RollingHash2D(matrix1)
  rh2 = RollingHash2D(matrix2)
  hash1 = rh1.get_hash(0,0,4,4)
  hash2 = rh2.get_hash(0,0,4,4)
  print(f"hash1: {hash1}")
  print(f"hash2: {hash2}")
  if hash1 == hash2:
      print("The two matrices are the same.")
  else:
      print("The two matrices are different.")
