function [systems] = generate_system(p)
systems(1) = tf([-1, 0], [1, p(1) + p(2), p(1) * p(2)]);
systems(2) = tf([-1], [1, p(1) + p(2), p(1) * p(2)]);
systems(3) = tf([1], [1, p(1) + p(2), p(1) * p(2)]);
end

