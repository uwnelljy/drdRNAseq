import random

# Set default parameters
TOTAL_NUMBER = 1000
GROUP_NUMBER = 2
LOWER = 0
UPPER = 1000


# Generate random relative expression level of two isoforms.
class Generate:
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def generate_iso_level(self):
        iso1 = random.randint(self.lower, self.upper)
        iso2 = random.randint(self.lower, self.upper)
        return iso1, iso2


class Main:
    def __init__(self):
        pass

    @staticmethod
    def write_txt():
        for i in range(TOTAL_NUMBER):
            for j in range(GROUP_NUMBER):
                iso1, iso2 = Generate(LOWER, UPPER).generate_iso_level()
                file = open('./expr{}_{}.txt'.format(i+1, j+1), 'w')
                file.write('iso1\t'+'{}'.format(iso1)+'\n')
                file.write('iso2\t'+'{}'.format(iso2))
                file.close()


if __name__ == '__main__':
    Main.write_txt()
