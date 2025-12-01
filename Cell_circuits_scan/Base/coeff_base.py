
import re

class coeff ():
    def __init__ (self, p):
        self.p = p
        self.check_is_format()
        self.symbol,self.letter,self.number1,self.number2 = self.split_string_1(self.p)
        self.number = str(self.number1)+str(self.number2)
        self.symbol,self.rest = self.split_symbol(self.p)
    def to_string(self):
        return self.p
    
    def check_is_format (self):
        pattern = r"([+-])([ab])(\d)(\d)"
        match = re.match(pattern, self.p)
        if not match:
            print("input doesnt match expected pattern")

    def split_string(self,s:str = None):
        if s==None:
            s = self.p
        pattern = r"([+-])([ab])(\d{2})"
        match = re.match(pattern, s)
        if match:
            symbol, letter, number = match.groups()
            return symbol, letter, number
        else:
            return None

    def split_string_1(self,s:str = None):
        if s==None:
            s = self.p
        pattern = r"([+-])([ab])(\d)(\d)"
        match = re.match(pattern, s)
        if match:
            symbol, letter, num1, num2 = match.groups()
            return symbol, letter, num1, num2
        else:
            return None

    def split_symbol(self,s:str = None):
        if s==None:
            s = self.p
        pattern = r"([+-])(.*)"
        match = re.match(pattern, s)
        if match:
            symbol, rest = match.groups()
            return symbol, rest
        else:
            return None

    