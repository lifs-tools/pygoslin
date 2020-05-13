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



class BaseParserEventHandler:

    def __init__(self):
        self.registered_events = {}
        self.rule_names = set()
        self.parser = None
        self.content = None
    
    
    
    # checking if all registered events are reasonable and orrur as rules in the grammar
    def sanity_check(self):
        
        for event_name in self.registered_events:
            if not event_name.endswith("_pre_event") and not event_name.endswith("_post_event"):
                raise Exception("Parser event handler error: event '%s' does not contain the suffix '_pre_event' or '_post_event'" % event_name);
            
            rule_name = event_name.replace("_pre_event", "").replace("_post_event", "")
            if rule_name not in self.rule_names:
                raise Exception("Parser event handler error: rule '%s' in event '%s' is not present in the grammar%s" % (rule_name, event_name, " '" + self.parser.grammar_name + "'" if self.parser != None else ""))
    
    
    def handle_event(self, event_name, node):
        if event_name in self.registered_events:
            self.registered_events[event_name](node)
            
