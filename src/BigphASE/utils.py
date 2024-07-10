import os
import time

__all__ = ['target_exist',
           'Record',
           'ExonTree',
           'make_dir',
           ]

#Record decorator 
def Record(func):
    def wrapper(*args, **kwargs):
        start = int(time.time())
        result = func(*args, **kwargs)
        end = int(time.time())
        setattr(result, "task", func.__name__)
        setattr(result, "start", start)
        setattr(result, "end", end)
        setattr(result, "duration", end - start)
        return result
    return wrapper

def make_dir(dir):
    if dir != '' and (not os.path.exists(dir)):
        os.makedirs(dir)
    else:
        return
    
def target_exist(*targets):
    """
    return True if all targets exist.
    """
    for target in targets:
        if  not os.path.exists(target):
            raise FileNotFoundError("File dos not exist:",target)
    return True


# ExonTree搜索树的节点
class ExonNode:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.height = 1
        self.left = None
        self.right = None

# 用于构建外显子区域的平衡二叉树，用于查找某个位点（代表位点的整数）是否落在外显子区域。
class ExonTree:
    def __init__(self):
        self.root = None

    def insert(self, start, end):
        node = ExonNode(start, end)
        if not self.root:
            self.root = node
            return

        stack = []
        current = self.root

        while current:
            stack.append(current)
            if start < current.start:
                if not current.left:
                    current.left = node
                    break
                current = current.left
            else:
                if not current.right:
                    current.right = node
                    break
                current = current.right

        while stack:
            current = stack.pop()
            current.height = 1 + max(self.get_height(current.left), self.get_height(current.right))
            balance = self.get_balance(current)

            if balance > 1:
                if start < current.left.start:
                    if stack:
                        parent = stack[-1]
                        if parent.left == current:
                            parent.left = self.right_rotate(current)
                        else:
                            parent.right = self.right_rotate(current)
                    else:
                        self.root = self.right_rotate(current)
                else:
                    current.left = self.left_rotate(current.left)
                    if stack:
                        parent = stack[-1]
                        if parent.left == current:
                            parent.left = self.right_rotate(current)
                        else:
                            parent.right = self.right_rotate(current)
                    else:
                        self.root = self.right_rotate(current)

            if balance < -1:
                if start > current.right.start:
                    if stack:
                        parent = stack[-1]
                        if parent.left == current:
                            parent.left = self.left_rotate(current)
                        else:
                            parent.right = self.left_rotate(current)
                    else:
                        self.root = self.left_rotate(current)
                else:
                    current.right = self.right_rotate(current.right)
                    if stack:
                        parent = stack[-1]
                        if parent.left == current:
                            parent.left = self.left_rotate(current)
                        else:
                            parent.right = self.left_rotate(current)
                    else:
                        self.root = self.left_rotate(current)

    def search(self, point) -> bool:
        current = self.root
        while current:
            if current.start <= point <= current.end:
                return True
            elif point < current.start:
                current = current.left
            else:
                current = current.right
        return False

    def get_height(self, node):
        if not node:
            return 0
        return node.height

    def get_balance(self, node):
        if not node:
            return 0
        return self.get_height(node.left) - self.get_height(node.right)

    def right_rotate(self, z):
        y = z.left
        if y is None:
            return z 

        T3 = y.right

        y.right = z
        z.left = T3

        z.height = 1 + max(self.get_height(z.left), self.get_height(z.right))
        y.height = 1 + max(self.get_height(y.left), self.get_height(y.right))

        return y

    def left_rotate(self, z):
        y = z.right
        if y is None:
            return z 

        T2 = y.left

        y.left = z
        z.right = T2

        z.height = 1 + max(self.get_height(z.left), self.get_height(z.right))
        y.height = 1 + max(self.get_height(y.left), self.get_height(y.right))

        return y