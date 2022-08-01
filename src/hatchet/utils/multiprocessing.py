import sys
import traceback
from multiprocessing import (
    Process,
    Queue,
    JoinableQueue,
    Lock,
    Value,
    cpu_count,
)
from hatchet.utils import ProgressBar as pb


class TaskHandler(Process):
    def __init__(self, handler_id, worker, task_queue, result_queue, progress_bar):
        Process.__init__(self)
        self.handler_id = handler_id
        self.worker = worker
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.progress_bar = progress_bar

    def run(self):
        while True:
            args = self.task_queue.get()
            if args is None:
                self.task_queue.task_done()
                break

            try:
                result = self.worker.work(*args)
            except Exception as e:
                typ, value, tb = sys.exc_info()
                result = RuntimeError(e).with_traceback(tb)
                # Traceback object is not preserved across multiprocessing; Save the traceback string as an attribute
                result.error = traceback.format_exception(typ, value, tb)

            if self.progress_bar is not None:
                self.progress_bar.progress(advance=True)

            self.result_queue.put((self.handler_id, result))
            self.task_queue.task_done()


class Worker:
    def __init__(self):
        pass

    def work(
        self,
        *args,
    ):
        raise NotImplementedError

    def run(self, work, n_instances=None, show_progress=True):

        n_work = len(work)
        if n_instances is None:
            n_instances = min(cpu_count(), n_work)
        else:
            n_instances = min(n_instances, n_work)

        if show_progress:
            progress_bar = pb.ProgressBar(
                total=n_work,
                length=40,
                lock=Lock(),
                counter=Value('i', 0),
                verbose=True,
            )
        else:
            progress_bar = None

        task_queue = JoinableQueue()
        for arg in work:
            if isinstance(arg, tuple):
                task_queue.put(arg)
            else:
                task_queue.put(tuple([arg]))

        for _ in range(n_instances):
            task_queue.put(None)

        result_queue = Queue()
        handlers = [TaskHandler(i, self, task_queue, result_queue, progress_bar) for i in range(n_instances)]

        for h in handlers:
            h.start()

        task_queue.join()

        try:
            _results = {}
            for handler_id, handler_result in [result_queue.get() for _ in range(n_work)]:
                if isinstance(handler_result, Exception):
                    error_string = ''.join(getattr(handler_result, 'error', []))
                    raise handler_result.__class__(f'HANDLER {handler_id} FAILED\n\n{error_string}')
                else:
                    _results[handler_id] = handler_result
            sorted_results = [_results[k] for k in sorted(_results.keys())]

        finally:
            task_queue.close()
            result_queue.close()

            for h in handlers:
                h.terminate()
                h.join()

        return sorted_results
