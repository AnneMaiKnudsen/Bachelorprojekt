from gwf import Workflow

gwf = Workflow()

gwf.target('MyTarget', inputs=[], outputs=['greeting.txt']) << """
echo hello world > greeting.txt
"""

gwf.target('TargetA', inputs=[], outputs=['x.txt']) << """
sleep 20 && echo "this is x" > x.txt
"""

gwf.target('TargetB', inputs=[], outputs=['y.txt']) << """
sleep 30 && echo "this is y" > y.txt
"""

gwf.target('TargetC', inputs=['x.txt', 'y.txt'], outputs=['z.txt']) << """
sleep 10 && cat x.txt y.txt > z.txt
"""
