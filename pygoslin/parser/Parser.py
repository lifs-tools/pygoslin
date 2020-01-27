from enum import Enum
from os import path

from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.GoslinFragmentParserEventHandler import GoslinFragmentParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from itertools import product as iter_product
from itertools import combinations as iter_combinations
import pygoslin


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
    
    
# DP stands for dynamic programming
class DPNode:
    def __init__(self, _rule1, _rule2):
        self.rule_index_1, self.left = _rule1
        self.rule_index_2, self.right = _rule2
        
        
        

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
        
        
        
        
     
def compute_rule_key(rule_index_1, rule_index_2):
    return (rule_index_1 << Parser.SHIFT) | rule_index_2
        
        
        
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
    
    
    
    
    
    def __init__(self, _parserEventHandler, grammar_filename, _quote = DEFAULT_QUOTE):
        self.next_free_rule_index = Parser.START_RULE
        self.TtoNT = {}
        self.NTtoNT = {}
        self.NTtoRule = {}
        self.originalNTtoNT = {}
        self.quote = _quote
        self.parser_event_handler = _parserEventHandler
        self.word_in_grammar = False
        self.grammar_name = ""
        self.used_eof = False
        
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
                        if c not in self.TtoNT: self.TtoNT[c] = set()
                        self.TtoNT[c].add(new_rule_index)
                    	
                    
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
        
        
        
        keys = set(self.TtoNT.keys())
        for c in keys:
            rules = set(self.TtoNT[c])
            self.TtoNT[c].clear()
            for rule in rules:
                for p in self.collect_one_backwards(rule):
                    
                    key = compute_rule_key(p, rule)
                    self.TtoNT[c].add(key)
        
        self.originalNTtoNT = {k: set(v) for k, v in self.NTtoNT.items()}
        
        keysNT = set(self.NTtoNT.keys())
        for r in keysNT:
            rules = set(self.NTtoNT[r])
            for rule in rules:
                for p in self.collect_one_backwards(rule): self.NTtoNT[r].add(p)
    
    
    
    
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
    def is_terminal(product_token, quote):
        return product_token[0] == quote and product_token[-1] == quote and len(product_token) > 2
    
    
    
    
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
            if c not in self.TtoNT: self.TtoNT[c] = set()
            next_index = self.get_next_free_rule_index()
            self.TtoNT[c].add(next_index)
            terminal_rules.append(next_index)
        
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
        if child_rule_index not in self.originalNTtoNT: return None
        
        
        for previous_index in self.originalNTtoNT[child_rule_index]:
            if previous_index == parent_rule_index:
                return []
            
            elif previous_index in self.originalNTtoNT:
                collection = self.collect_backwards(previous_index, parent_rule_index)
                if collection != None:
                    collection.append(previous_index)
                    return collection
        return None
        
        
        
        
    
    def raise_events(self, node):
        if node != None:
            node_rule_name = self.NTtoRule[node.rule_index] if node.fire_event else ""
            if node.fire_event: self.parser_event_handler.handle_event(node_rule_name + "_pre_event", node)
            
            if node.left != None: # node.terminal is != None when node is leaf
                self.raise_events(node.left)
                if node.right != None: self.raise_events(node.right)
                
            if node.fire_event: self.parser_event_handler.handle_event(node_rule_name + "_post_event", node)
            
    
    
    
    

    # filling the syntax tree including events
    def fill_tree(self, node, dp_node):
        # checking and extending nodes for single rule chains
        key = compute_rule_key(dp_node.rule_index_1, dp_node.rule_index_2) if dp_node.left != None else dp_node.rule_index_2
        
        merged_rules = self.collect_backwards(key, node.rule_index)
        if merged_rules != None:
            for rule_index in merged_rules:
                node.left = TreeNode(rule_index, rule_index in self.NTtoRule)
                node = node.left
        
        
        if dp_node.left != None: # None => leaf
            node.left = TreeNode(dp_node.rule_index_1, dp_node.rule_index_1 in self.NTtoRule)
            node.right = TreeNode(dp_node.rule_index_2, dp_node.rule_index_2 in self.NTtoRule)
            self.fill_tree(node.left, dp_node.left)
            self.fill_tree(node.right, dp_node.right)
            
        else:
            # I know, it is not 100% clean to store the character in an integer
            # especially when it is not the dedicated attribute for, but the heck with it!
            node.terminal = dp_node.rule_index_1
    
    
    
    
    # re-implementation of Cocke-Younger-Kasami algorithm
    def parse(self, text_to_parse):
        if self.used_eof: text_to_parse += Parser.EOF_SIGN
        
        self.parse_regular(text_to_parse)
        return self.parser_event_handler.content
        
        
        
    def parse_regular(self, text_to_parse):
        self.word_in_grammar = False
        self.parser_event_handler.content = None
        
        n = len(text_to_parse)
        # dp stands for dynamic programming, nothing else
        dp_table = [None for x in range(n)]
        # Ks is a lookup, which fields in the dp_table are filled
        Ks = [set() for x in range(n)]

        
        
        
        
        
        for i in range(n):
            dp_table[i] = [None for x in range(n - i)]
            for j in range(n - i): dp_table[i][j] = {}
        
        for i in range(n):
            c = text_to_parse[i]
            if c not in self.TtoNT: return
            
            for rule_index in sorted(self.TtoNT[c]):
                new_key = rule_index >> Parser.SHIFT
                old_key = rule_index & Parser.MASK
                #dp_node = DPNode(c, old_key, None, None)
                dp_node = DPNode([c, None], [old_key, None])
                dp_table[i][0][new_key] = dp_node
                Ks[i].add(0)
                
                
        sft = Parser.SHIFT
        nt = self.NTtoNT
        
        for i in range (1, n):
            im1 = i - 1
            for j in range(n - i):
                D, jp1 = dp_table[j], j + 1
                Di = D[i]
                
                adding = False
                
                for k in Ks[j]:
                    jpok = jp1 + k
                    D1, D2 = D[k], dp_table[jpok][im1 - k]
                    if im1 - k not in Ks[jpok]: continue
                
                    for index_pair_1, index_pair_2 in iter_product(D1, D2):
                        key = (index_pair_1 << sft) | index_pair_2
                        
                        if key not in nt: continue
                        
                        content = DPNode([index_pair_1, D1[index_pair_1]], [index_pair_2, D2[index_pair_2]])
                        for rule_index in nt[key]: Di[rule_index] = content
                
                if len(D[i]) > 0: Ks[j].add(i)
        
        for i in range(n - 1, 0, -1):
            if Parser.START_RULE in dp_table[0][i]:
                self.word_in_grammar = True
                parse_tree = TreeNode(Parser.START_RULE, Parser.START_RULE in self.NTtoRule)
                self.fill_tree(parse_tree, dp_table[0][i][Parser.START_RULE])
                self.raise_events(parse_tree)
                break
        
    
    
class GoslinParser(Parser):
    def __init__(self):
        self.event_handler = GoslinParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/Goslin.g4", Parser.DEFAULT_QUOTE)
        
        
        
class GoslinFragmentParser(Parser):
    def __init__(self):
        self.event_handler = GoslinFragmentParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/GoslinFragments.g4", Parser.DEFAULT_QUOTE)
        
        
class LipidMapsParser(Parser):
    def __init__(self):
        self.event_handler = LipidMapsParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/LipidMaps.g4", Parser.DEFAULT_QUOTE)
        
        
class LipidParser:
    def __init__(self):
        self.parser_list = [GoslinParser(), GoslinFragmentParser(), LipidMapsParser()]
        
    def parse(self, lipid_name):
        
        for parser in self.parser_list:
            lipid = parser.parse(lipid_name)
            if lipid != None:
                return lipid
            
        return None