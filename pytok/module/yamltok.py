import yaml
from yaml import *
import os.path

'''
New loader definition from
	https://stackoverflow.com/questions/528281/how-can-i-include-an-yaml-file-inside-another#15437697
which allows to include yaml files in yaml files
'''
class LoaderTok(yaml.Loader):
	def __init__(self, stream, additonalAnchors = None):
		self._file_root_dir = os.path.split(stream.name)[0]
		self._found_anchors = {}
		super(LoaderTok, self).__init__(stream)
		if additonalAnchors is not None:
			self.anchors.update(additonalAnchors)


	def compose_node(self, parent, index):
		ret = yaml.Loader.compose_node(self, parent, index)
		self._found_anchors.update(self.anchors)
		return ret;

	def include(self, node):
		filename = os.path.join(self._file_root_dir , self.construct_scalar(node))
		with open(filename, 'r') as f:
			return yaml.load(f, lambda s: LoaderTok(s,self._found_anchors))

LoaderTok.add_constructor('!include', LoaderTok.include)


# override resolver to avoid on/off/yes/no is converted to a boolean
# from https://stackoverflow.com/questions/36463531/pyyaml-automatically-converting-certain-keys-to-boolean-values
for ch in "OoYyNn":
	if len(yaml.resolver.Resolver.yaml_implicit_resolvers[ch]) == 1:
		del yaml.resolver.Resolver.yaml_implicit_resolvers[ch]
	else:
		yaml.resolver.Resolver.yaml_implicit_resolvers[ch] = [
		                   x for x in yaml.resolver.Resolver.yaml_implicit_resolvers[ch]
		                     if x[0] != 'tag:yaml.org,2002:bool']


####################################################3
# merge results
# inspired and partly copied from https://gist.github.com/aryzhov/d92991528cc7a996711a220b37903756
def _trailing_plus_count(s):
	for i in range(0, len(s)):
		if s[-i-1] != '+':
			return i
	return len(s)


def _get_base_key(key):
	return key.rstrip("+")


def _merge_list_or_dict(v, v_extend):
	if type(v) is list and type(v_extend) is list:
		v.extend(v_extend)
	elif type(v) is dict and type(v_extend) is dict:
		v.update(v_extend)
	else:
		raise TypeError("Can merge only lists or dicts")
	return v


def merge_lists(root):
	if type(root) == dict:
		to_merge = set()
		for k in root.keys():
			if type(k) == str and _trailing_plus_count(k) > 0:
				to_merge.add(_get_base_key(k))
		for base_key in to_merge:
			keys_sorted = sorted( [ k for k in root.keys()
			                          if type(k) == str and _get_base_key(k) == base_key ],
			                      key = _trailing_plus_count )

			v = None
			for k in keys_sorted:
				if v is None:
					v = root[k]
				else:
					_merge_list_or_dict(v, root[k])
				root.pop(k)
			root[base_key] = v

		for v in root.values():
			if type(v) is list or type(v) is dict:
				merge_lists(v)
	elif type(root) is list:
		for node in root:
			merge_lists(node)
