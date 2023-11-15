from collections import namedtuple
import os
from pathlib import Path
import tempfile
from unittest import TestCase
from unittest.mock import MagicMock
from . import guanin


class TestOutput:
    def test_default_output(self, monkeypatch):
        monkeypatch.setattr(os, "name", "posix")
        app_dir = guanin.default_output(make_dir=False)

        assert app_dir == Path(os.getenv("XDG_DATA_HOME")) / "guanin"

    def test_create_dirs(self, monkeypatch):
        mock = MagicMock(return_value=True)
        monkeypatch.setattr(Path, "mkdir", mock)

        guanin.default_output(make_dir=True)

        assert mock.called_with(parents=True, exist_ok=True)


class TestSetup(TestCase):
    def setUp(self):
        self.Args = namedtuple("Args", ["output"])
        self.output_dir = Path(__file__).parent / "test_files"

        self.original_posix_env = os.environ.get("XDG_DATA_HOME")
        os.environ["XDG_DATA_HOME"] = str(Path(__file__).parent / "test_files")
        self.original_nt_env = os.environ.get("LOCALAPPDATA")
        os.environ["LOCALAPPDATA"] = str(Path(__file__).parent / "test_files")

    def tearDown(self):
        if self.original_posix_env:
            os.environ["XDG_DATA_HOME"] = self.original_posix_env
        if self.original_nt_env:
            os.environ["LOCALAPPDATA"] = self.original_nt_env
        try:
            Path(self.output_dir / "guanin").rmdir()
        except FileNotFoundError:
            pass

    def test_defined_output(self):
        mock = MagicMock(return_value="posix")
        os.name = mock()

        args = self.Args(output=self.output_dir)
        run = guanin.Setup(args)

        assert run.output == self.output_dir
        assert run.output.exists()

    def test_default_output_posix(self):
        mock = MagicMock(return_value="posix")
        os.name = mock()

        args = self.Args(output="")
        run = guanin.Setup(args)

        assert run.output == self.output_dir / "guanin"

    def test_default_output_nt(self):
        mock = MagicMock(return_value="nt")
        os.name = mock

        args = self.Args(output="")
        run = guanin.Setup(args)

        assert run.output == self.output_dir / "guanin"

    def test_output_is_created(self):
        mock = MagicMock(return_value="posix")
        os.name = mock

        args = self.Args(output="")
        run = guanin.Setup(args)

        assert Path(run.output).exists()

    def test_default_output_posix_no_environment(self):
        mock = MagicMock(return_value="posix")
        os.name = mock()

        os.environ.pop("XDG_DATA_HOME")

        args = self.Args(output="")
        run = guanin.Setup(args)

        assert run.output == Path(tempfile.gettempdir()) / os.getlogin() \
            / "guanin"
