from datetime import datetime
from math import floor
from time import sleep, time


class ProcessBar:

    def __init__(self, n: int, width: int = 50) -> None:
        self.time = time()
        self.n = n
        self.width = width
        return None

    def plot(self, x: int) -> None:
        deltaTime = time() - self.time
        x = x / self.n
        deltaTime *= (1 - x) / (x + 1e-3)
        hours, deltaTime = divmod(deltaTime, 3600)
        minutes, seconds = divmod(deltaTime, 60)
        bar = f'|{"█" * floor(x * self.width):<{self.width}}|'
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> {bar} ({hours:02.0f}:{minutes:02.0f}:{seconds:02.0f}).', end = '\r' if x < 1 else '\n', flush = True)
        return None


if __name__ == '__main__':
    processBar = ProcessBar(100)
    for i in range(1, 101):
        sleep(0.1)
        processBar.plot(i)
