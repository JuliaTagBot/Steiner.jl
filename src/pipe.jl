using Sockets, Random
import JSON2

function handle_pipe_input!(sock, input)
   try
      output = handle_solve_conics(JSON2.read(input))
      write(sock, JSON2.write(output), "\n")
   catch err
      @warn err
   end
end


"""
   sync_setup_pipe(pipe_name)

Setup a named pipe to communicate with the conic solver. This will run async.
"""
function setup_pipe(pipe_name::String)
    @async begin
      server = listen(pipe_name)
       while true
          sock = accept(server)
          @async while isopen(sock)
             handle_pipe_input!(sock, readline(sock))
          end
       end
    end
end

"""
   sync_setup_pipe(pipe_name)

Setup a named pipe to communicate with the conic solver. This will block the current session.
"""
function sync_setup_pipe(pipe_name::String)
    @sync begin
      @async wait(setup_pipe(pipe_name))
    end
end
