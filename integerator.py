class Integrator:
    def __init__(self, f, h):
        self.f = f
        self.h = h

    def euler(self):
        x = [self.x0]
        y = [self.y0]

        while True:
            next_y = y[-1] + self.h * self.f(x[-1], y[-1])
            y.append(next_y)
            x.append(x[-1] + self.h)
            if abs(y[-1] - y[-2]) > 1e100:
                break

        return x, y
