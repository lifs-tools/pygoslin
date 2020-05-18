"""
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
                   Nils Hoffmann  -  nils.hoffmann {at} isas.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


from enum import Enum
from os import path

import pygoslin
from pygoslin.domain.LipidExceptions import LipidParsingException, LipidException
from itertools import combinations as iter_combinations
from pygoslin.parser.SumFormulaParserEventHandler import SumFormulaParserEventHandler


pyx_support = True
try:
    import pyximport
    from pygoslin.parser.ParserCore import parser_core
except:
    pyx_support = False


class Context(Enum):
    NoContext = 1
    InLineComment = 2
    InLongComment = 3
    InQuote = 4
    
class MatchWords(Enum):
    NoMatch = 1
    LineCommentStart = 2
    LineCommentEnd = 3
    LongCommentStart = 4
    LongCommentEnd = 5
    Quote = 6
    
    
        

class TreeNode:

    def __init__(self, _rule, _fire_event):
        self.rule_index = _rule
        self.left = None
        self.right = None
        self.terminal = 0
        self.fire_event = _fire_event
    
    def get_text(self):
        if self.terminal == 0:
            left_str = self.left.get_text()
            right_str = self.right.get_text() if self.right != None else ""
            return "%s%s" % (left_str if left_str != Parser.EOF_SIGN else "", right_str if right_str != Parser.EOF_SIGN else "")
        return self.terminal
        
        
        
        
     
compute_rule_key = lambda rule_index_1, rule_index_2: (rule_index_1 << Parser.SHIFT) | rule_index_2
        
        
        
class Parser:
    
    SHIFT = 32
    MASK = (1 << SHIFT) - 1
    RULE_ASSIGNMENT = ':'
    RULE_SEPARATOR = '|'
    RULE_TERMINAL = ';'
    EOF_RULE_NAME = "EOF"
    EOF_SIGN = chr(1)
    EOF_RULE = 1
    START_RULE = 2
    DEFAULT_QUOTE = "'"
    
    
    
    
    def get_next_free_rule_index(self):
        if self.next_free_rule_index <= Parser.MASK:
            n = self.next_free_rule_index
            self.next_free_rule_index += 1
            return n
        raise Exception("Error: grammar is too big.")
    
    
    
    
    
    def __init__(self, _parserEventHandler, grammar_filename = None, _quote = DEFAULT_QUOTE):
        self.next_free_rule_index = Parser.START_RULE
        self.TtoNT = {}
        self.NTtoNT = {}
        self.NTtoRule = {}
        self.quote = _quote
        self.parser_event_handler = _parserEventHandler
        self.word_in_grammar = False
        self.grammar_name = ""
        self.used_eof = False
        self.substitution = {}
        self.OTtoNT = {}
        self.right_pair = []
        self.left_pair = []
        
        if grammar_filename == None: return
        
        if path.exists(grammar_filename):
            # interpret the rules and create the structure for parsing
            rules = Parser.extract_text_based_rules(grammar_filename, self.quote)
            self.grammar_name = Parser.split_string(rules[0], ' ', self.quote)[1]
            del rules[0]
            ruleToNT = {Parser.EOF_RULE_NAME: Parser.EOF_RULE}
            self.TtoNT[Parser.EOF_SIGN] = set([Parser.EOF_RULE])
            for rule_line in rules:
                
                tokens_level_1 = []
                for t in Parser.split_string(rule_line, Parser.RULE_ASSIGNMENT, self.quote): tokens_level_1.append(t.strip(' '))
                if len(tokens_level_1) != 2: raise Exception("Error: corrupted token in grammar rule: '%s'" % rule_line);
                
                if len(Parser.split_string(tokens_level_1[0], ' ', self.quote)) > 1:
                    raise Exception("Error: several rule names on left hand side in grammar rule: '%s'" % rule_line);
                

                rule = tokens_level_1[0]
                
                if rule ==  Parser.EOF_RULE_NAME:
                    raise Exception("Error: rule name is not allowed to be called EOF")
                
                products = [p.strip(' ') for p in Parser.split_string(tokens_level_1[1], Parser.RULE_SEPARATOR, self.quote)]
                
                if rule not in ruleToNT: ruleToNT[rule] = self.get_next_free_rule_index()
                new_rule_index = ruleToNT[rule]
                
                if new_rule_index not in self.NTtoRule: self.NTtoRule[new_rule_index] = rule
                
                
                for product in products:
                    non_terminals, non_terminal_rules = [], []
                    for NT in Parser.split_string(product, ' ', self.quote):
                        stripedNT = NT.strip(' ')
                        if Parser.is_terminal(stripedNT, self.quote): stripedNT = Parser.de_escape(stripedNT, self.quote)
                        non_terminals.append(stripedNT)
                        self.used_eof |= stripedNT == Parser.EOF_RULE_NAME
                        
                    
                    NTFirst = non_terminals[0]
                    if len(non_terminals) > 1 or not Parser.is_terminal(NTFirst, self.quote) or len(NTFirst) != 3:
                        for non_terminal in non_terminals:
                            
                            if Parser.is_terminal(non_terminal, self.quote):
                                non_terminal_rules.append(self.add_terminal(non_terminal))
                                
                            else:
                                if non_terminal not in ruleToNT:
                                    ruleToNT[non_terminal] = self.get_next_free_rule_index()
                                non_terminal_rules.append(ruleToNT[non_terminal])
                                
                    else:
                        c = NTFirst[1]
                        t_rule = 0
                        if c not in self.TtoNT: 
                            t_rule = self.get_next_free_rule_index()
                            self.TtoNT[c] = set([t_rule])
                        else:
                            t_rule = list(self.TtoNT[c])[0]
                            
                        if t_rule not in self.NTtoNT: self.NTtoNT[t_rule] = set()
                        self.NTtoNT[t_rule].add(new_rule_index)
                    	
                    
                    # more than two rules, insert intermediate rule indexes
                    while len(non_terminal_rules) > 2:
                        rule_index_2 = non_terminal_rules.pop()
                        rule_index_1 = non_terminal_rules.pop()
                        
                        key = compute_rule_key(rule_index_1, rule_index_2)
                        next_index = self.get_next_free_rule_index()
                        if key not in self.NTtoNT: self.NTtoNT[key] = set()
                        self.NTtoNT[key].add(next_index)
                        non_terminal_rules.append(next_index)
                    
                        
                    # two product rules
                    if len(non_terminal_rules) == 2:
                        rule_index_2 = non_terminal_rules[1]
                        rule_index_1 = non_terminal_rules[0]
                        key = compute_rule_key(rule_index_1, rule_index_2)
                        if key not in self.NTtoNT: self.NTtoNT[key] = set()
                        self.NTtoNT[key].add(new_rule_index)
                    
                    # only one product rule
                    elif len(non_terminal_rules) == 1:
                        rule_index_1 = non_terminal_rules[0]
                        if rule_index_1 == new_rule_index:
                            raise Exception("Error: corrupted token in grammar: rule '%s' is not allowed to refer soleley to itself." % rule)
                        
                        if rule_index_1 not in self.NTtoNT: self.NTtoNT[rule_index_1] = set()
                        self.NTtoNT[rule_index_1].add(new_rule_index)
            
            # adding all rule names into the event handler
            for rule_name in ruleToNT:
                self.parser_event_handler.rule_names.add(rule_name)
                
            self.parser_event_handler.parser = self
            self.parser_event_handler.sanity_check()
            
        else:
            raise Exception("Error: file '%s' does not exist or can not be opened." % grammar_filename)
        
        def top_nodes( rule_index):
            collection, collection_top = [rule_index], []
            i = 0
            while i < len(collection):
                current_index = collection[i]
                if current_index in self.NTtoNT:
                    for previous_index in self.NTtoNT[current_index]: collection.append(previous_index)
                else:
                    collection_top.append(current_index)
                i += 1
            return collection_top
        
        
        
        
        for rule_index, values in self.NTtoNT.items():
            for rule in values:
                for rule_top in top_nodes(rule):
                    chain = self.collect_backwards(rule, rule_top)
                    if chain == None: continue
                    chain = [rule_top] + chain + [rule]
                    while len(chain) > 1:
                        top = chain[0]
                        chain = chain[1:]
                        self.substitution[rule_index + (top << 16)] = chain
                    
        for k, v in self.TtoNT.items():
            self.OTtoNT[k] = list(v)[0]
        
        
        for d in [self.TtoNT, self.NTtoNT]:
            for rules in d.values():
                new_rules = set()
                for rule in rules:
                    new_rules |= set([p for p in self.collect_one_backwards(rule)])
                rules |= new_rules
                    
            
    
        self.right_pair = [set() for i in range(self.next_free_rule_index)]
        self.left_pair = [set() for i in range(self.next_free_rule_index)]
        for key in self.NTtoNT:
            if key <= self.MASK: continue
            self.right_pair[key >> self.SHIFT].add(key)
            self.left_pair[key & self.MASK].add(key)
           
           
           
    
    def extract_text_based_rules(grammar_filename, quote = DEFAULT_QUOTE):
        grammar, sb,current_position = "", [], 0
        
        with open(grammar_filename, mode = "rt") as infile:
            grammar = infile.read() + "\n";
        grammar_length = len(grammar)
        
        # deleting comments to prepare for splitting the grammar in rules.
        # Therefore, we have to consider three different contexts, namely
        # within a quote, within a line comment, within a long comment.
        # As long as we are in one context, key words for starting / ending
        # the other contexts have to be ignored.
        current_context = Context.NoContext
        last_escaped_backslash = -1
        for i in range (grammar_length - 1):
            match = MatchWords.NoMatch
            
            if i > 0 and grammar[i] == '\\' and grammar[i - 1] == '\\' and last_escaped_backslash != i - 1:
                last_escaped_backslash = i
                continue
            
            if grammar[i] == '/' and grammar[i + 1] == '/': match = MatchWords.LineCommentStart
            elif grammar[i] == '\n': match = MatchWords.LineCommentEnd
            elif grammar[i] == '/' and grammar[i + 1] == '*': match = MatchWords.LongCommentStart
            elif grammar[i] == '*' and grammar[i + 1] == '/': match = MatchWords.LongCommentEnd
            elif grammar[i] == quote and not (i >= 1 and grammar[i - 1] == '\\' and i - 1 != last_escaped_backslash): match = MatchWords.Quote
            
            if match != MatchWords.NoMatch:
                
                if current_context == Context.NoContext:
                    if match == MatchWords.LongCommentStart:
                        sb.append(grammar[current_position : i])
                        current_context = Context.InLongComment
                        
                    elif match == MatchWords.LineCommentStart:
                        sb.append(grammar[current_position : i])
                        current_context = Context.InLineComment
                        
                    elif match == MatchWords.Quote:
                        current_context = Context.InQuote
                    
                elif current_context == Context.InQuote:
                    if match == MatchWords.Quote:
                        current_context = Context.NoContext;
                    
                    
                elif current_context == Context.InLineComment:
                    if match == MatchWords.LineCommentEnd:
                        current_context = Context.NoContext
                        current_position = i + 1
                    
                elif current_context == Context.InLongComment:
                    if match == MatchWords.LongCommentEnd:
                        current_context = Context.NoContext
                        current_position = i + 2
                
                    
                    
        if current_context == Context.NoContext:
            sb.append(grammar[current_position : grammar_length])
            
        else:
            raise Exception("Error: corrupted grammar '%s', ends either in comment or quote" % grammar_filename)
        
        grammar = "".join(sb).replace("\r\n", "").replace("\n", "").replace("\r", "").strip(" ");
        if grammar[-1] != Parser.RULE_TERMINAL:
            raise Exception("Error: corrupted grammar'%s', last rule has no termininating sign, was: '%s'" % (grammar_filename, grammar[-1]))
        
        rules = Parser.split_string(grammar, Parser.RULE_TERMINAL, quote)
        
        if len(rules) < 1:
            raise Exception("Error: corrupted grammar '%s', grammar is empty" % grammar_filename)
        
        grammar_name_rule = Parser.split_string(rules[0], ' ', quote)
        
        if grammar_name_rule[0] != "grammar":
            raise Exception("Error: first rule must start with the keyword 'grammar'")
        
        elif len(grammar_name_rule) != 2:
            raise Exception("Error: incorrect first rule")
        
        return rules
    
    
    
    
    def dump_child(self, class_name, parser_event_handler_name):
        with open("%s.py" % class_name, mode = "wt") as dumpfile:
            dumpfile.write("from pygoslin.parser.Parser import Parser\n")
            dumpfile.write("from pygoslin.parser.%s import %s\n\n" % (parser_event_handler_name, parser_event_handler_name))
            dumpfile.write("class %s(Parser):\n" % class_name)
            dumpfile.write("    def __init__(self):\n")
            dumpfile.write("        super().__init__(%s(), None, Parser.DEFAULT_QUOTE)\n" % parser_event_handler_name)
            dumpfile.write("        self.next_free_rule_index = %s\n" % str(self.next_free_rule_index))
            dumpfile.write("        self.TtoNT = %s\n" % str(self.TtoNT))
            dumpfile.write("        self.NTtoNT = %s\n" % str(self.NTtoNT))
            dumpfile.write("        self.NTtoRule = %s\n" % str(self.NTtoRule))
            dumpfile.write("        self.word_in_grammar = %s\n" % str(self.word_in_grammar))
            dumpfile.write("        self.used_eof = %s\n" % str(self.used_eof))
            dumpfile.write("        self.substitution = %s\n" % str(self.substitution))
            dumpfile.write("        self.OTtoNT = %s\n" % str(self.OTtoNT))
            dumpfile.write("        self.right_pair = %s\n" % str(self.right_pair))
            dumpfile.write("        self.left_pair = %s\n" % str(self.left_pair))
    
    
    
    def split_string(text, separator, quote = DEFAULT_QUOTE):
        in_quote, tokens, sb, last_char, last_escaped_backslash = False, [], [], '\0', False
        
        for c in text:
            escaped_backslash = False
            if not in_quote:
                if c == separator:
                    if len(sb) > 0: tokens.append("".join(sb))
                    sb = []
                else:
                    if c == quote: in_quote = not in_quote
                    sb.append(c)
                    
            else:
                if c == '\\' and last_char == '\\' and not last_escaped_backslash:
                    escaped_backslash = True
                    
                elif c == quote and not (last_char == '\\' and not last_escaped_backslash):
                    in_quote = not in_quote
                sb.append(c)
            last_escaped_backslash = escaped_backslash
            last_char = c
                
        if len(sb) > 0:
            tokens.append("".join(sb))
        if in_quote: raise Exception("Error: corrupted token in grammar")
        
        return tokens
    
    
    
    # checking if string is terminal
    is_terminal = lambda product_token, quote: product_token[0] == quote and product_token[-1] == quote and len(product_token) > 2
    
    
    
    
    def de_escape(text, quote):
        # remove the escape chars
        sb, last_escape_char = [], False
        for c in text:
            escape_char = False
            
            if c != '\\':
                sb.append(c)
                
            else:
                if not last_escape_char: escape_char = True
                else: sb.append(c)
            
            last_escape_char = escape_char
        
        return "".join(sb)
    
    
    
    # splitting the whole terminal in a tree structure where characters of terminal are the leafs and the inner nodes are added non terminal rules
    def add_terminal(self, text):
        terminal_rules = []
        for i in range(1, len(text) - 1):
            c = text[i]
            t_rule = 0
            if c not in self.TtoNT:
                t_rule = self.get_next_free_rule_index()
                self.TtoNT[c] = set([t_rule])
            else:
                t_rule = list(self.TtoNT[c])[0]
            terminal_rules.append(t_rule)
        
        while len(terminal_rules) > 1:
            rule_index_2 = terminal_rules.pop()
            rule_index_1 = terminal_rules.pop()
            
            next_index = self.get_next_free_rule_index()
            
            key = compute_rule_key(rule_index_1, rule_index_2)
            if key not in self.NTtoNT: self.NTtoNT[key] = set()
            self.NTtoNT[key].add(next_index)
            terminal_rules.append(next_index)
        
        return terminal_rules[0]
    
    
    
    
    # expanding singleton rules, e.g. S -> A, A -> B, B -> C
    def collect_one_backwards(self, rule_index):
        collection = [rule_index]
        i = 0
        while i < len(collection):
            current_index = collection[i]
            if current_index in self.NTtoNT:
                for previous_index in self.NTtoNT[current_index]: collection.append(previous_index)
            i += 1
        return collection
    
    
    
    
    def collect_backwards(self, child_rule_index, parent_rule_index):
        if child_rule_index not in self.NTtoNT: return None
        
        for previous_index in self.NTtoNT[child_rule_index]:
            if previous_index == parent_rule_index:
                return []
            
            elif previous_index in self.NTtoNT:
                collection = self.collect_backwards(previous_index, parent_rule_index)
                if collection != None:
                    collection.append(previous_index)
                    return collection
        return None
        
        
        
        
    
    def raise_events(self, node):
        if node:
            node_rule_name = self.NTtoRule[node.rule_index] if node.fire_event else ""
            if node.fire_event: self.parser_event_handler.handle_event(node_rule_name + "_pre_event", node)
            
            if node.left: # node.terminal is != None when node is leaf
                self.raise_events(node.left)
                if node.right: self.raise_events(node.right)
                
            if node.fire_event: self.parser_event_handler.handle_event(node_rule_name + "_post_event", node)
            
    
    
    
    

    # filling the syntax tree including events
    def fill_tree(self, node, parse_content):
        # checking and extending nodes for single rule chains
        
        if parse_content[1] != None:
            bottom_rule = compute_rule_key(parse_content[0], parse_content[2])
            top_rule = node.rule_index
        else:
            top_rule = parse_content[2]
            bottom_rule = self.OTtoNT[parse_content[0]]
        
        
        s_key = bottom_rule + (top_rule << 16)
        if bottom_rule != top_rule and s_key in self.substitution:
            for rule_index in self.substitution[s_key]:
                node.left = TreeNode(rule_index, rule_index in self.NTtoRule)
                node = node.left
        
        
        if parse_content[1] != None: # None => leaf
            node.left = TreeNode(parse_content[0], parse_content[0] in self.NTtoRule)
            node.right = TreeNode(parse_content[2], parse_content[2] in self.NTtoRule)
            self.fill_tree(node.left, parse_content[1])
            self.fill_tree(node.right, parse_content[3])
            
        else:
            node.terminal = parse_content[0]
    
    
    
    
    def parse(self, text_to_parse, raise_error = True):
        old_lipid = text_to_parse
        if self.used_eof: text_to_parse += Parser.EOF_SIGN
        
        self.parse_regular(text_to_parse)
        if raise_error and not self.word_in_grammar:
            raise LipidParsingException("Lipid '%s' can not be parsed by grammar '%s'" % (old_lipid, self.grammar_name))
        
        return self.parser_event_handler.content
        
        
        
        
        
        
    def parse_regular(self, text_to_parse):
        self.word_in_grammar = False
        self.parser_event_handler.content = None
        
        # call parser core function written in cython for performance boost
        dp_table = parser_core(text_to_parse, self.NTtoNT, self.TtoNT, self.left_pair, self.right_pair) if pyx_support else self.parser_pure(text_to_parse)
        if dp_table == None: return
        
        
        for i in range(len(text_to_parse) - 1, 0, -1):
            if Parser.START_RULE in dp_table[i][0]:
                self.word_in_grammar = True
                parse_tree = TreeNode(Parser.START_RULE, Parser.START_RULE in self.NTtoRule)
                self.fill_tree(parse_tree, dp_table[i][0][Parser.START_RULE])
                self.raise_events(parse_tree)
                break
            if self.used_eof: break
            
            
            
            
            
    # re-implementation of Cocke-Younger-Kasami algorithm
    def parser_pure(self, text_to_parse):
        n = len(text_to_parse)
        rgt = self.right_pair
        lft = self.left_pair
        
        
        DP = [[{} for j in range(n - i)] for i in range(n)]
        DL = [[set() for j in range(n - i)] for i in range(n)]
        DR = [[set() for j in range(n - i)] for i in range(n)]
        Ks = [set([0]) for i in range(n)]
        
        t, nt, mask, shift = self.TtoNT, self.NTtoNT, (1 << 32) - 1, 32
            
        for i in range(n):
            c = text_to_parse[i]
            if c not in t: return None
            
            for rule_index in t[c]:
                DP[0][i][rule_index] = [c, None, rule_index, None]
                DL[0][i] |= lft[rule_index]
                DR[0][i] |= rgt[rule_index]
        
        
        for i in range (1, n):
            im1, DPi = i - 1, DP[i]
            for j in range(n - i):
                jp1, DPij, Ksj = j + 1, DPi[j], Ks[j]
                
                for k in Ksj:
                    jpok, im1mk = jp1 + k, im1 - k
                    if im1mk not in Ks[jpok]: continue
                    
                    for key in DR[k][j] & DL[im1mk][jpok]:
                        index_pair_1 = key >> shift
                        index_pair_2 = key & mask
                        parse_content = [index_pair_1, DP[k][j][index_pair_1], index_pair_2, DP[im1mk][jpok][index_pair_2]]
                        for rule_index in nt[key]:
                            DPij[rule_index] = parse_content
                            DL[i][j] |= lft[rule_index]
                            DR[i][j] |= rgt[rule_index]
                                    
                if DPij: Ksj.add(i)

        return DP
           
           
class SumFormulaParser(Parser):
    def __init__(self):
        self.event_handler = SumFormulaParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/SumFormula.g4", Parser.DEFAULT_QUOTE)

        
