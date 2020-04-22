function discNode = MountOnNode(Position,nodePos)
for i=1:length(Position)
    discNode(i) = find(min(abs(Position(i)-nodePos))==abs(Position(i)-nodePos));
end
end
