def xgcd(x, y):
    old_r, r = (x, y)
    old_s, s = (1, 0)
    old_t, t = (0, 1)

    while r != 0:
        quotient = old_r // r
        old_r, r = (r, old_r - quotient * r)
        old_s, s = (s, old_s - quotient * s)
        old_t, t = (t, old_t - quotient * t)

    return old_s, old_t, old_r  # a, b, g


class FieldElement:
    def __init__(self, value):
        field = Field()
        self.value = value % field.p
        self.field = field

    def __add__(self, right):
        return self.field.add(self, right)

    def __mul__(self, right):
        return self.field.multiply(self, right)

    def __sub__(self, right):
        return self.field.subtract(self, right)

    def __truediv__(self, right):
        return self.field.divide(self, right)

    def __neg__(self):
        return self.field.negate(self)

    def inverse(self):
        return self.field.inverse(self)

    # modular exponentiation -- be sure to encapsulate in parentheses!
    def __xor__(self, exponent):
        acc = FieldElement(1)
        val = FieldElement(self.value)
        for i in reversed(range(len(bin(exponent)[2:]))):
            acc = acc * acc
            if (1 << i) & exponent != 0:
                acc = acc * val
        return acc

    def __eq__(self, other):
        if isinstance(other, FieldElement):
            return self.value == other.value
        else:
            return False

    def __neq__(self, other):
        return self.value != other.value

    def __repr__(self):
        return str(self.value)

    def __bytes__(self):
        return bytes(str(self).encode())

    def __index__(self):
        return self.value

    def __hex__(self):
        return hex(self.value)

    def is_zero(self):
        if self.value == 0:
            return True
        else:
            return False


class Field:
    def __init__(self):
        self.p = 2**64 - 2**32 + 1

    def zero(self):
        return FieldElement(0)

    def one(self):
        return FieldElement(1)

    # MULTIPLICATIVE_GROUP_GENERATOR
    def g(self):
        return FieldElement(7)

    def multiply(self, left, right):
        return FieldElement((left.value * right.value) % self.p)

    def add(self, left, right):
        return FieldElement((left.value + right.value) % self.p)

    def subtract(self, left, right):
        return FieldElement((self.p + left.value - right.value) % self.p)

    def negate(self, operand):
        return FieldElement((self.p - operand.value) % self.p)

    def inverse(self, operand):
        a, b, g = xgcd(operand.value, self.p)
        return FieldElement(a)

    def divide(self, left, right):
        assert (not right.is_zero()), "divide by zero"
        a, b, g = xgcd(right.value, self.p)
        return FieldElement(left.value * a % self.p)
