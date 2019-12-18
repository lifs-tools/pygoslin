
class BaseParserEventHandler:

    def __init__(self):
        self.registered_events = {}
        self.rule_names = set()
        self.parser = None
    
    
    
    # checking if all registered events are reasonable and orrur as rules in the grammar
    def sanity_check(self):
        
        for event_name in self.registered_events:
            if not event_name.endswith("_pre_event") and not event_name.endswith("_post_event"):
                raise Exception("Parser event handler error: event '%s' does not contain the suffix '_pre_event' or '_post_event'" % event_name);
            
            rule_name = event_name.replace("_pre_event", "").replace("_post_event", "")
            if rule_name not in self.rule_names:
                raise Exception("Parser event handler error: rule '%s' in event '%s' is not present in the grammar%s" + (rule_name, event_name, " '" + self.parser.grammar_name + "'" if self.parser != None else ""))
    
    
    def handle_event(self, event_name, node):
        if event_name in self.registered_events:
            self.registered_events[event_name](node)
            