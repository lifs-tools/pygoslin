
class BaseParserEventHandler:

    def __init__(self):
        self.registered_events = []
        rule_names = set()
        parser = None
    
    
    
    # checking if all registered events are reasonable and orrur as rules in the grammar
    def sanityCheck(self):
        
        for event_name in registered_events:
            if not event_name.endswith("_pre_event") and event_name.endswith("_post_event")):
                raise Exception("Parser event handler error: event '%s' does not contain the suffix '_pre_event' or '_post_event'" % event_name);
            
            rule_name = event_name.replace("_pre_event", "").replace("_post_event", "")
            if rule_name not in rule_names:
                throw new Exception("Parser event handler error: rule '%s' in event '%s' is not present in the grammar" + (rule_name, event_name  " '" + parser.grammarName + "'" if parser != None else ""))
    
    
    del handleEvent(event_name, node):
        if event_name in registeredEvents:
            registeredEvents[eventName](node)
            