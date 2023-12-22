import pgscatalog.corelib
import pgscatalog.calclib

def run():
    test = pgscatalog.corelib.ScoringFile()
    test2 = pgscatalog.calclib.TestClass()
    print(f"I imported a test object {test}")
    print(f"I imported a test object {test2} from a different package")
    
if __name__ == "__main__":
    run()