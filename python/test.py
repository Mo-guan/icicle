class myInt:
    def __init__(self, value):
        self.value = value

    def __add__(self, other):
        return myInt(self.value + other.value)

    def __repr__(self):
        return "value = " + str(self.value)


if __name__ == "__main__":
    a = [myInt(1), myInt(2), myInt(3)]
    b = a[0:2]
    c = a[0:1]
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")

    b[0].value = 10
    print(f"a = {a}")

    c[0].value = 20
    print(f"a = {a}")
